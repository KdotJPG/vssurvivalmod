﻿using System.Collections.Generic;
using Vintagestory.API.Client;
using Vintagestory.API.Common;
using Vintagestory.API.Datastructures;
using Vintagestory.API.MathTools;

namespace Vintagestory.GameContent
{
    public class BlockBehaviorWrenchOrientable : BlockBehavior
    {
        public static Dictionary<string, List<AssetLocation>> VariantsByType = new Dictionary<string, List<AssetLocation>>();

        public string BaseCode;
        bool hideInteractionHelpInSurvival;
        private static List<ItemStack> wrenchItems = new List<ItemStack>();

        public BlockBehaviorWrenchOrientable(Block block) : base(block)
        {
        }

        public override WorldInteraction[] GetPlacedBlockInteractionHelp(IWorldAccessor world, BlockSelection selection, IPlayer forPlayer, ref EnumHandling handling)
        {
            if (hideInteractionHelpInSurvival && forPlayer?.WorldData.CurrentGameMode == EnumGameMode.Survival) return base.GetPlacedBlockInteractionHelp(world, selection, forPlayer, ref handling);

            handling = EnumHandling.PassThrough;
            if (wrenchItems.Count == 0)   // This is a potentially rather slow wildcard search of all items (especially if mods add many items) therefore we want to run this only once per game
            {
                Item[] wrenches = world.SearchItems(new AssetLocation("wrench-*"));
                foreach(Item item in wrenches) wrenchItems.Add(new ItemStack(item));
            }
            if (wrenchItems.Count > 0)
            {
                return new WorldInteraction[] { new WorldInteraction()
                {
                    ActionLangCode = "Rotate",
                    Itemstacks = wrenchItems.ToArray(),
                    MouseButton = EnumMouseButton.Right
                } };
            }
            else return new WorldInteraction[0];
        }

        public override void Initialize(JsonObject properties)
        {
            base.Initialize(properties);

            hideInteractionHelpInSurvival = properties["hideInteractionHelpInSurvival"].AsBool(false);
            BaseCode = properties["baseCode"].AsString();
        }

        public override void OnLoaded(ICoreAPI api)
        {
            base.OnLoaded(api);

            if (BaseCode != null)
            {
                List<AssetLocation> vars;
                if (!VariantsByType.TryGetValue(BaseCode, out vars)) VariantsByType[BaseCode] = vars = new List<AssetLocation>();
                vars.Add(block.Code);
            }
        }

        public override void OnUnloaded(ICoreAPI api)
        {
            wrenchItems.Clear();
        }
    }
}
