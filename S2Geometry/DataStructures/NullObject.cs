using System;
using System.ComponentModel;

namespace Google.Common.Geometry
{
    internal struct NullObject<T> : IEquatable<T>
    {
        [DefaultValue(false)]
        private bool isNotNull;// default property initializers are not supported for structs

        private NullObject(T item, bool isnull) : this()
        {
            this.isNotNull = isNotNull;
            this.Item = item;
        }

        public NullObject(T item) : this(item, item == null)
        {
        }

        public static NullObject<T> Null()
        {
            return new NullObject<T>();
        }

        public T Item { get; private set; }

        public bool IsNull()
        {
            return !this.isNotNull;
        }

        public bool IsNotNull()
        {
            return this.isNotNull;
        }

        public static implicit operator T(NullObject<T> nullObject)
        {
            return nullObject.Item;
        }

        public static implicit operator NullObject<T>(T item)
        {
            return new NullObject<T>(item);
        }

        public override string ToString()
        {
            return (Item != null) ? Item.ToString() : "NULL";
        }

        public override int GetHashCode()
        {
            unchecked
            {
                if (this.IsNull())
                return 0;

                var result = Item.GetHashCode();

                if (result >= 0)
                    result++;
            
                return result;
            }
        }

        public bool Equals(T obj)
        {
            if (obj == null)
                return this.IsNull();

            if (!(obj is T))
                return false;

            var no = (NullObject<T>)obj;

            if (this.IsNull())
                return no.IsNull();

            if (no.IsNull())
                return false;

            return this.Item.Equals(no.Item);
        }
    }
}
