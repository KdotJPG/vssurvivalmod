﻿using Vintagestory.API;
using Vintagestory.API.Common;
using Vintagestory.API.Common.Entities;
using Vintagestory.API.MathTools;

namespace Vintagestory.GameContent
{
    public class BlockBehaviorHeatSource : BlockBehavior, IHeatSource
    {
        float heatStrength;

        public BlockBehaviorHeatSource(Block block) : base(block)
        {
        }


        public override void Initialize(JsonObject properties)
        {
            heatStrength = properties["heatStrength"].AsFloat(0);

            base.Initialize(properties);
        }

        public float GetHeatStrength(IWorldAccessor world, BlockPos heatSourcePos, BlockPos heatReceiverPos)
        {
            if (block.EntityClass != null)
            {
                var behs = world.BlockAccessor.GetBlockEntity(heatSourcePos) as IHeatSource;
                if (behs != null)
                {
                    return behs.GetHeatStrength(world, heatSourcePos, heatReceiverPos) * GameMath.Clamp(1 - (heatSourcePos.DistanceTo(heatReceiverPos) - 1)/3f, 0, 1);
                }
            }

            return heatStrength * GameMath.Clamp(1 - (heatSourcePos.DistanceTo(heatReceiverPos) - 1) / 3f, 0, 1);
        }
    }
}
