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

    sealed class KeyValuePairEqualityComparer<K, V> : IEqualityComparer<KeyValuePair<K, V>>
    {
        readonly IEqualityComparer<K> keyequalityComparer;


        /// <summary>
        /// Create an entry equalityComparer using the default equalityComparer for keys
        /// </summary>
        public KeyValuePairEqualityComparer() { keyequalityComparer = EqualityComparer<K>.Default; }


        /// <summary>
        /// Create an entry equalityComparer from a specified item equalityComparer for the keys
        /// </summary>
        /// <param name="keyequalityComparer">The key equalityComparer</param>
        public KeyValuePairEqualityComparer(IEqualityComparer<K> keyequalityComparer)
        {
            if (keyequalityComparer == null)
                throw new NullReferenceException("Key equality comparer cannot be null");
            this.keyequalityComparer = keyequalityComparer;
        }


        /// <summary>
        /// Get the hash code of the entry
        /// </summary>
        /// <param name="entry">The entry</param>
        /// <returns>The hash code of the key</returns>
     
        public int GetHashCode(KeyValuePair<K, V> entry) { return keyequalityComparer.GetHashCode(entry.Key); }


        /// <summary>
        /// Test two entries for equality
        /// </summary>
        /// <param name="entry1">First entry</param>
        /// <param name="entry2">Second entry</param>
        /// <returns>True if keys are equal</returns>
        public bool Equals(KeyValuePair<K, V> entry1, KeyValuePair<K, V> entry2)
        { return keyequalityComparer.Equals(entry1.Key, entry2.Key); }
    }
}
