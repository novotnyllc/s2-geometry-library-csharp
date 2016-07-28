using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry.DataStructures
{
    // Adapted from https://github.com/sestoft/C5/blob/master/C5/hashing/HashBag.cs
    class HashBag<T> : ICollection<T>
    {
        readonly Dictionary<T, int> dict;
        int size;

        public HashBag() : this(EqualityComparer<T>.Default)
        {
        }

        public HashBag(IEqualityComparer<T> itemEqualityComparer)
        {
            dict = new Dictionary<T, int>(itemEqualityComparer);
        }

        public IEnumerator<T> GetEnumerator()
        {
            foreach (var item in dict)
            {
                for (var i = 0; i < item.Value; i++)
                    yield return item.Key;
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public void Add(T item)
        {
            int val;
            if (dict.TryGetValue(item, out val))
            {
                dict[item] = ++val;
            }
            else
            {
                dict.Add(item, 1);
            }
            size++;
        }

        public void Clear()
        {
            dict.Clear();
            size = 0;
        }

        public bool Contains(T item)
        {
            return dict.ContainsKey(item);
        }

        public void CopyTo(T[] array, int arrayIndex)
        {
            if (arrayIndex < 0 || arrayIndex + Count > array.Length)
                throw new ArgumentOutOfRangeException();

            foreach (var p in dict)
                for (var j = 0; j < p.Value; j++)
                    array[arrayIndex++] = p.Key;
        }

        public bool Remove(T item)
        {
            int val;

            if (dict.TryGetValue(item, out val))
            {
                size--;
                if (val == 1)
                    dict.Remove(item);
                else
                {
                    dict[item] = --val;
                }

                return true;
            }
            return false;
        }

        public int Count => size;
        public bool IsReadOnly => false;
    }
}
