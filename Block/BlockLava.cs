using System;
using System.Collections.Generic;
using Vintagestory.API.Client;
using Vintagestory.API.Common;
using Vintagestory.API.Common.Entities;
using Vintagestory.API.MathTools;

namespace Vintagestory.GameContent
{
    /// <summary>
    /// Handles lava burning nearby combustible blocks. Searches y levels above the lava in a radius
    /// of 2,3,2 respectively. If the block is combustible and has a burn temperature lower than or equal
    /// to the temperature at that location then a fire block will be placed in the adjacent air block.
    /// </summary>
    public class BlockLava : BlockForFluidsLayer, IBlockFlowing
    {
        public string Flow { get; set; }
        public Vec3i FlowNormali { get => null; set {} }
        public bool IsLava => true;

        /// <summary>
        /// Data structure returned to the tick system to be used by this block in order to
        /// initialize the BEFire with the right BlockFacing value.
        /// </summary>
        private class FireLocation
        {
            public readonly BlockPos firePos;
            public readonly BlockFacing facing;

            public FireLocation(BlockPos firePos, BlockFacing facing)
            {
                this.firePos = firePos;
                this.facing = facing;
            }
        }
        /// <summary>
        /// Temperature of lava. Controls determining whether an item should burn
        /// </summary>
        private readonly int temperature = 1200;

        /// <summary>
        /// Amount of temperature is decreased for each one block distance away from lava(Manhattan distance)
        /// </summary>
        private readonly int tempLossPerMeter = 100;

        private Block blockFire;
        AdvancedParticleProperties[] fireParticles;

        public BlockLava() : base()
        {
            if (Attributes != null)
            {
                temperature = Attributes["temperature"].AsInt(1200);
                tempLossPerMeter = Attributes["tempLossPerMeter"].AsInt(100);
            }
        }

        public override void OnLoaded(ICoreAPI api)
        {
            base.OnLoaded(api);

            if (blockFire == null)
            {
                blockFire = api.World.GetBlock(new AssetLocation("fire"));

                fireParticles = new AdvancedParticleProperties[blockFire.ParticleProperties.Length];
                for (int i = 0; i < fireParticles.Length; i++)
                {
                    fireParticles[i] = blockFire.ParticleProperties[i].Clone();
                }

                fireParticles[2].HsvaColor[2].avg += 60;
                fireParticles[2].LifeLength.avg += 3;
            }
        }

        public override void OnServerGameTick(IWorldAccessor world, BlockPos pos, object extra = null)
        {
            base.OnServerGameTick(world, pos, extra);

            FireLocation fireLocation = (FireLocation)extra;
            world.BlockAccessor.SetBlock(blockFire.BlockId,fireLocation.firePos);
            BlockEntity befire = world.BlockAccessor.GetBlockEntity(fireLocation.firePos);
            befire?.GetBehavior<BEBehaviorBurning>()?.OnFirePlaced(fireLocation.facing, null);
        }

        /// <summary>
        /// Searches for an air block next to a combustible block
        /// </summary>
        /// <param name="world"></param>
        /// <param name="lavaPos"></param>
        /// <returns>The position of the air block next to a combustible block</returns>
        private FireLocation SearchAreaForAirNextToCombustibleBlock(IWorldAccessor world, BlockPos lavaPos)
        {
            FireLocation combustibleBlockPos = SearchRadiusForAirNextToCombustibleBlock(world, lavaPos, 1, 2);
            if(combustibleBlockPos == null)
            {
                combustibleBlockPos = SearchRadiusForAirNextToCombustibleBlock(world, lavaPos, 2, 3);
                if(combustibleBlockPos == null)
                {
                    combustibleBlockPos = SearchRadiusForAirNextToCombustibleBlock(world, lavaPos, 3, 2);
                }
            }
            return combustibleBlockPos;
        }

        /// <summary>
        /// Searches a given horizontal radius for an air block next to a combustible block
        /// </summary>
        /// <param name="world"></param>
        /// <param name="lavaPos"></param>
        /// <param name="y">Current y level</param>
        /// <param name="radius">Horizontal Radius</param>
        /// <returns></returns>
        private FireLocation SearchRadiusForAirNextToCombustibleBlock(IWorldAccessor world, BlockPos lavaPos, int y, int radius)
        {
            for (int x = -radius; x <= radius; x++)
            {
                for (int z = -radius; z <= radius; z++)
                {
                    BlockPos airBlockPos = lavaPos.AddCopy(x, y, z);
                    
                    BlockFacing facing = IsNextToCombustibleBlock(world, lavaPos, airBlockPos);
                    if (facing != null)
                    {
                        return new FireLocation(airBlockPos, facing);
                    }
                }
            }
            return null;
        }

        /// <summary>
        /// Returns true if the given air block position is next to a combustible block, false otherwise. The
        /// block must be combustible and at a burnable temperature
        /// </summary>
        /// <param name="world"></param>
        /// <param name="lavaPos"></param>
        /// <param name="airBlockPos"></param>
        /// <returns></returns>
        private BlockFacing IsNextToCombustibleBlock(IWorldAccessor world, BlockPos lavaPos, BlockPos airBlockPos)
        {
            Block airBlock = world.BlockAccessor.GetBlock(airBlockPos);
            if (airBlock.BlockId == 0)
            {
                foreach (BlockFacing facing in BlockFacing.ALLFACES)
                {
                    if (facing != BlockFacing.DOWN)
                    {
                        BlockPos combustibleBlockPos = airBlockPos.Copy().Add(facing);
                        Block combustibleBlock = world.BlockAccessor.GetBlock(combustibleBlockPos);
                        if (combustibleBlock.CombustibleProps != null && combustibleBlock.CombustibleProps.BurnTemperature <= GetTemperatureAtLocation(lavaPos, airBlockPos))
                        {
                            return facing;
                        }
                    }
                }
            }
            return null;
        }

        /// <summary>
        /// Returns the temperature at the given location based on it's distance from the lava.
        /// </summary>
        /// <param name="lavaPos"></param>
        /// <param name="airBlockPos"></param>
        /// <returns></returns>
        private int GetTemperatureAtLocation(BlockPos lavaPos, BlockPos airBlockPos)
        {
            int distance = lavaPos.ManhattenDistance(airBlockPos);
            return temperature - (distance * tempLossPerMeter);
        }

        public override bool ShouldReceiveServerGameTicks(IWorldAccessor world, BlockPos pos, Random offThreadRandom, out object extra)
        {
            if (LiquidLevel == 7)
            {                
                FireLocation fireLocation = SearchAreaForAirNextToCombustibleBlock(world, pos);
                if(fireLocation != null)
                {
                    extra = fireLocation;
                    return true;
                }
            }

            extra = null;
            return false;
        }




        public override bool ShouldReceiveClientParticleTicks(IWorldAccessor world, IPlayer player, BlockPos pos, out bool isWindAffected)
        {
            isWindAffected = false;
            Block block = world.BlockAccessor.GetBlock(pos.X, pos.Y + 1, pos.Z);
            return !block.IsLiquid() && (block.CollisionBoxes == null || block.CollisionBoxes.Length == 0);
        }

        public override void OnAsyncClientParticleTick(IAsyncParticleManager manager, BlockPos pos, float windAffectednessAtPos, float secondsTicking)
        {
            if (GameMath.MurmurHash3Mod(pos.X, pos.Y, pos.Z, 100) < 2)
            {
                for (int i = 0; i < fireParticles.Length; i++)
                {
                    AdvancedParticleProperties bps = fireParticles[i];
                    bps.Quantity.avg = i * 0.3f; // No cubes, medium fire, double smoke
                    bps.WindAffectednesAtPos = windAffectednessAtPos;
                    bps.basePos.X = pos.X + TopMiddlePos.X;
                    bps.basePos.Y = pos.Y + TopMiddlePos.Y;
                    bps.basePos.Z = pos.Z + TopMiddlePos.Z;

                    manager.Spawn(bps);
                }
            }

            base.OnAsyncClientParticleTick(manager, pos, windAffectednessAtPos, secondsTicking);
        }
    }
}
