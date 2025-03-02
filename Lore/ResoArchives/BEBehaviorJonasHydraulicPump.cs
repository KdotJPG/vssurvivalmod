﻿using System.Text;
using Vintagestory.API.Client;
using Vintagestory.API.Common;
using Vintagestory.API.Config;
using Vintagestory.API.Datastructures;
using Vintagestory.API.Util;

namespace Vintagestory.GameContent
{
    public class BEBehaviorJonasHydraulicPump : BEBehaviorControlPointAnimatable, INetworkedLight
    {
        ControlPoint cp;
        bool on;
        string networkCode;
        public string oncommands;
        public string offcommands;

        bool hasTempGear;
        bool hasFlywheel;

        bool IsRepaired => hasTempGear && hasFlywheel;
        bool ReceivesPower => (animControlPoint?.ControlData as AnimationMetaData)?.AnimationSpeed > 0;

        protected override Shape AnimationShape {
            get {
                AssetLocation shapePath = Block.ShapeInventory.Base.Clone().WithPathPrefixOnce("shapes/").WithPathAppendixOnce(".json");
                var shape = Shape.TryGet(Api, shapePath);
                return shape;
            }
        }

        public BEBehaviorJonasHydraulicPump(BlockEntity blockentity) : base(blockentity)
        {
        }


        public override void Initialize(ICoreAPI api, JsonObject properties)
        {
            this.Api = api;
            registerToLightNetworkControlPoint();
            base.Initialize(api, properties);
            updatePumpingState();
        }

        public void setNetwork(string networkCode)
        {
            this.networkCode = networkCode;
            registerToLightNetworkControlPoint();
            Blockentity.MarkDirty(true);
            updatePumpingState();
        }

        void registerToLightNetworkControlPoint()
        {
            if (networkCode == null) return;

            modSys = Api.ModLoader.GetModSystem<ModSystemControlPoints>();
            var controlpointcode = AssetLocation.Create(networkCode, Block.Code.Domain);
            cp = modSys[controlpointcode];

            on = IsRepaired && ReceivesPower;
        }

        protected override void BEBehaviorControlPointAnimatable_Activate(ControlPoint cpoint)
        {
            if (IsRepaired)
            {
                base.BEBehaviorControlPointAnimatable_Activate(cpoint);
            } else
            {
                this.animControlPoint = cpoint;
            }

            on = IsRepaired && ReceivesPower;
            updatePumpingState();
        }

        internal void Interact(IPlayer byPlayer, BlockSelection blockSel)
        {
            if (Api.Side != EnumAppSide.Server) return;

            if (byPlayer.WorldData.CurrentGameMode == EnumGameMode.Creative && byPlayer.InventoryManager.ActiveHotbarSlot.Empty)
            {
                if (byPlayer.Entity.Controls.Sneak)
                {
                    hasTempGear = false;
                    hasFlywheel = false;
                } else
                {
                    on = !on;
                    updatePumpingState();
                }

                Blockentity.MarkDirty(true);
                return;
            }

            if (!IsRepaired)
            {
                var slot = byPlayer.InventoryManager.ActiveHotbarSlot;
                if (slot.Itemstack?.Collectible.Code.Path == "largegear-temporal" && !hasTempGear)
                {
                    if (byPlayer.WorldData.CurrentGameMode != EnumGameMode.Creative) slot.TakeOut(1);
                    hasTempGear = true;
                    Api.World.PlaySoundAt(new AssetLocation("sounds/effect/latch"), Pos.X + 0.5, Pos.Y, Pos.Z + 0.5, null, true, 16);
                }
                if (slot.Itemstack?.Collectible.Code.Path == "jonasparts-flywheel" && !hasFlywheel)
                {
                    if (byPlayer.WorldData.CurrentGameMode != EnumGameMode.Creative) slot.TakeOut(1);
                    hasFlywheel = true;
                    Api.World.PlaySoundAt(new AssetLocation("sounds/effect/latch"), Pos.X + 0.5, Pos.Y, Pos.Z + 0.5, null, true, 16);
                }

                on = IsRepaired && ReceivesPower;
                updatePumpingState();

                Blockentity.MarkDirty(true);

                return;
            }
        }

        protected void updatePumpingState()
        {
            if (Api == null) return;

            if (IsRepaired) updateAnimationstate();

            if (cp != null)
            {
                cp.ControlData = on;
                cp.Trigger();
            }

            var caller = new Caller()
            {
                CallerPrivileges = new string[] { "*" },
                Pos = Pos.ToVec3d(),
                Type = EnumCallerType.Block
            };

            var becs = Blockentity as BlockEntityCommands;
            if (becs != null)
            {
                becs.CallingPrivileges = caller.CallerPrivileges;
                becs.Execute(caller, on ? oncommands : offcommands);
            }            

            Blockentity.MarkDirty(true);
        }

        public override void FromTreeAttributes(ITreeAttribute tree, IWorldAccessor worldAccessForResolve)
        {
            base.FromTreeAttributes(tree, worldAccessForResolve);

            hasFlywheel = tree.GetBool("hasFlywheel");
            hasTempGear = tree.GetBool("hasTempGear");
            on = tree.GetBool("on");
            networkCode = tree.GetString("networkCode");
            if (networkCode == "") networkCode = null;
            oncommands = tree.GetString("oncommands");
            offcommands = tree.GetString("offcommands");

            updatePumpingState();
        }

        public override void ToTreeAttributes(ITreeAttribute tree)
        {
            base.ToTreeAttributes(tree);

            tree.SetBool("hasTempGear", hasTempGear);
            tree.SetBool("hasFlywheel", hasFlywheel);
            tree.SetBool("on", on);
            tree.SetString("networkCode", networkCode);
            tree.SetString("oncommands", oncommands);
            tree.SetString("offcommands", offcommands);
        }

        public override void GetBlockInfo(IPlayer forPlayer, StringBuilder dsc)
        {
            base.GetBlockInfo(forPlayer, dsc);

            if (Api is ICoreClientAPI capi && capi.Settings.Bool["extendedDebugInfo"] == true)
            {
                dsc.AppendLine("network code: " + networkCode);
                dsc.AppendLine(on ? "On" : "Off");
                dsc.AppendLine("oncommand:" + oncommands);
                dsc.AppendLine("offcommand:" + offcommands);
            }

            if (!ReceivesPower) dsc.AppendLine(Lang.Get("No power."));
            if (!hasTempGear) dsc.AppendLine(Lang.Get("Missing large temporal gear."));
            if (!hasFlywheel) dsc.AppendLine(Lang.Get("Missing pump head."));
        }


        public override bool OnTesselation(ITerrainMeshPool mesher, ITesselatorAPI tessThreadTesselator)
        {
            if (animUtil.activeAnimationsByAnimCode.Count == 0)
            {
                if (hasFlywheel)
                {
                    mesher.AddMeshData(genMesh(new AssetLocation("shapes/block/machine/jonas/pumphead1-flywheel.json")));
                }
                if (hasTempGear)
                {
                    mesher.AddMeshData(genMesh(new AssetLocation("shapes/block/machine/jonas/pumphead1-gear.json")));
                }

                return false;
            }
            else
            {
                return base.OnTesselation(mesher, tessThreadTesselator);
            }
        }

        private MeshData genMesh(AssetLocation assetLocation)
        {
            return ObjectCacheUtil.GetOrCreate(Api, "hydrpumpmesh-" + assetLocation.Path + "-" + Block.Shape.rotateY, () =>
            {
                var shape = Api.Assets.TryGet(assetLocation).ToObject<Shape>();
                (Api as ICoreClientAPI).Tesselator.TesselateShape(Block, shape, out var mesh, Block.Shape.RotateXYZCopy);
                return mesh;
            });
        }
    }
}
