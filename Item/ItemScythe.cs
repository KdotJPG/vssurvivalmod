﻿using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Vintagestory.API;
using Vintagestory.API.Client;
using Vintagestory.API.Common;
using Vintagestory.API.Common.Entities;
using Vintagestory.API.Config;
using Vintagestory.API.Datastructures;
using Vintagestory.API.MathTools;
using Vintagestory.API.Server;
using Vintagestory.API.Util;

namespace Vintagestory.GameContent
{
    public class ItemScythe : ItemShears
    {
        string[] allowedPrefixes;
        string[] disallowedSuffixes;
        SkillItem[] skillItems;

        public override int MultiBreakQuantity { get { return 5; } }

        public override void OnLoaded(ICoreAPI api)
        {
            base.OnLoaded(api);
            allowedPrefixes = Attributes["codePrefixes"].AsArray<string>();
            disallowedSuffixes = Attributes["disallowedSuffixes"].AsArray<string>();

            skillItems = new SkillItem[]
            {
                new SkillItem()
                {
                    Code = new AssetLocation("trim grass"),
                    Name = Lang.Get("Trim grass")
                }.WithLetterIcon(api as ICoreClientAPI, "T"),
                new SkillItem()
                {
                    Code = new AssetLocation("remove grass"),
                    Name = Lang.Get("Remove grass")
                }.WithLetterIcon(api as ICoreClientAPI, "R"),
            };
        }

        
        public override WorldInteraction[] GetHeldInteractionHelp(ItemSlot inSlot)
        {
            return base.GetHeldInteractionHelp(inSlot).Append(new WorldInteraction()
            {
                ActionLangCode = "heldhelp-settoolmode",
                HotKeyCode = "toolmodeselect"
            });
        }

        public override bool CanMultiBreak(Block block)
        {
            for (int i = 0; i < allowedPrefixes.Length; i++)
            {
                if (block.Code.Path.StartsWith(allowedPrefixes[i]))
                {
                    // Disable scything on thick snow variants (-snow2, -snow3 etc.)
                    if (disallowedSuffixes != null)
                    {
                        for (int j = 0; j < disallowedSuffixes.Length; j++)
                        {
                            if (block.Code.Path.EndsWith(disallowedSuffixes[j])) return false;
                        }
                    }

                    return true;
                }
            }

            return false;   
        }

        public override void OnHeldAttackStart(ItemSlot slot, EntityAgent byEntity, BlockSelection blockSel, EntitySelection entitySel, ref EnumHandHandling handling)
        {
            base.OnHeldAttackStart(slot, byEntity, blockSel, entitySel, ref handling);
            if (blockSel == null) return;

            IPlayer byPlayer = (byEntity as EntityPlayer)?.Player;
            if (!byEntity.World.Claims.TryAccess(byPlayer, blockSel.Position, EnumBlockAccessFlags.BuildOrBreak))
            {
                return;
            }

            byEntity.Attributes.SetBool("didBreakBlocks", false);
            byEntity.Attributes.SetBool("didPlayScytheSound", false);
            handling = EnumHandHandling.PreventDefault;
        }

        public override bool OnHeldAttackStep(float secondsPassed, ItemSlot slot, EntityAgent byEntity, BlockSelection blockSelection, EntitySelection entitySel)
        {
            if (byEntity.World.Side == EnumAppSide.Client)
            {
                ModelTransform tf = new ModelTransform();
                tf.EnsureDefaultValues();

                float t = secondsPassed / 1.35f;

                float f = (float)Easings.EaseOutBack(Math.Min(t * 2f, 1));
                float f2 = (float)Math.Sin(GameMath.Clamp(Math.PI * 1.4f * (t - 0.5f), 0, 3));

                tf.Translation.X += Math.Min(0.2f, t * 3);
                tf.Translation.Y -= Math.Min(0.75f, t * 3);
                tf.Translation.Z -= Math.Min(1, t * 3);
                tf.ScaleXYZ += Math.Min(1, t * 3);
                tf.Origin.X -= Math.Min(0.75f, t * 3);
                tf.Rotation.X = -Math.Min(30, t * 30) + f * 30 + (float)f2 * 120f;
                tf.Rotation.Z = -f * 110;

                if (secondsPassed > 1.75f)
                {
                    float b = 2 * (secondsPassed - 1.75f);
                    tf.Rotation.Z += b * 140;
                    tf.Rotation.X /= (1 + b * 10);
                    tf.Translation.X -= b * 0.4f;
                    tf.Translation.Y += b * 2 / 0.75f;
                    tf.Translation.Z += b * 2;
                }

                byEntity.Controls.UsingHeldItemTransformBefore = tf;
            }

            performActions(secondsPassed, byEntity, slot, blockSelection);

            // Crappy fix to make harvesting not buggy T_T
            if (api.Side == EnumAppSide.Server) return true;

            return secondsPassed < 2f;
        }

        public override void OnHeldAttackStop(float secondsPassed, ItemSlot slot, EntityAgent byEntity, BlockSelection blockSelection, EntitySelection entitySel)
        {
            performActions(secondsPassed, byEntity, slot, blockSelection);
        }


        bool trimMode;

        private void performActions(float secondsPassed, EntityAgent byEntity, ItemSlot slot, BlockSelection blockSelection)
        {
            if (blockSelection == null) return;

            var block = api.World.BlockAccessor.GetBlock(blockSelection.Position);
            IPlayer byPlayer = (byEntity as EntityPlayer)?.Player;

            var canmultibreak = CanMultiBreak(api.World.BlockAccessor.GetBlock(blockSelection.Position));

            if (canmultibreak && secondsPassed > 0.75f && byEntity.Attributes.GetBool("didPlayScytheSound") == false)
            {
                api.World.PlaySoundAt(new AssetLocation("sounds/tool/scythe1"), byEntity, byPlayer, true, 16);
                byEntity.Attributes.SetBool("didPlayScytheSound", true);
            }

            if (canmultibreak && secondsPassed > 1.05f && byEntity.Attributes.GetBool("didBreakBlocks") == false)
            {
                if (byEntity.World.Side == EnumAppSide.Server && byEntity.World.Claims.TryAccess(byPlayer, blockSelection.Position, EnumBlockAccessFlags.BuildOrBreak))
                {
                    trimMode = block.Variant["tallgrass"] != null && block.Variant["tallgrass"] != "eaten" && slot.Itemstack.Attributes.GetInt("toolMode", 0) == 0;

                    OnBlockBrokenWith(byEntity.World, byEntity, slot, blockSelection);
                }

                byEntity.Attributes.SetBool("didBreakBlocks", true);
            }
        }


        protected override void breakMultiBlock(BlockPos pos, IPlayer plr)
        {
            if (trimMode)
            {
                var block = api.World.BlockAccessor.GetBlock(pos);
                var trimmedBlock = api.World.GetBlock(block.CodeWithVariant("tallgrass", "eaten"));
                if (block == trimmedBlock) return;

                if (trimmedBlock != null)
                {
                    api.World.BlockAccessor.BreakBlock(pos, plr);
                    api.World.BlockAccessor.MarkBlockDirty(pos);

                    api.World.BlockAccessor.SetBlock(trimmedBlock.BlockId, pos);
                    var be = api.World.BlockAccessor.GetBlockEntity(pos) as BlockEntityTransient;
                    if (be != null) be.ConvertToOverride = block.Code.ToShortString();

                    return;
                }
            }

            base.breakMultiBlock(pos, plr);
        }


        public override SkillItem[] GetToolModes(ItemSlot slot, IClientPlayer forPlayer, BlockSelection blockSel)
        {
            return skillItems;
        }

        public override void SetToolMode(ItemSlot slot, IPlayer byPlayer, BlockSelection blockSelection, int toolMode)
        {
            slot.Itemstack.Attributes.SetInt("toolMode", toolMode);
        }

        public override int GetToolMode(ItemSlot slot, IPlayer byPlayer, BlockSelection blockSelection)
        {
            return slot.Itemstack.Attributes.GetInt("toolMode", 0);
        }


        public override void OnUnloaded(ICoreAPI api)
        {
            for (int i = 0; skillItems != null && i < skillItems.Length; i++)
            {
                skillItems[i]?.Dispose();
            }
        }


    }
}
