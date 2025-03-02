﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Vintagestory.API.Common;
using Vintagestory.API.MathTools;

namespace Vintagestory.GameContent
{
    public class TradeItem : JsonItemStack
    {
        //public string Name;
        public NatFloat Price;
        public NatFloat Stock;
        public RestockOpts Restock = new RestockOpts()
        {
            HourDelay = 24,
            Quantity = 1
        };
        public SupplyDemandOpts SupplyDemand = new SupplyDemandOpts()
        {
            PriceChangePerDay = 0.1f,
            PriceChangePerPurchase = 0.1f
        };

        public ResolvedTradeItem Resolve(IWorldAccessor world)
        {
            this.Resolve(world, "TradeItem");

            return new ResolvedTradeItem()
            {
                Stack = this.ResolvedItemstack,
                Price = (int)Math.Max(1, Math.Round(Price.nextFloat(1f, world.Rand))),
                Stock = Stock == null ? 0 : (int)Math.Round(Stock.nextFloat(1f, world.Rand)),
                Restock = Restock,
                SupplyDemand = SupplyDemand
            };
        }
    }
}
