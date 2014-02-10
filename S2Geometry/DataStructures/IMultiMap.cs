using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace System.Collections.Generic
{
    internal interface IMultiMap<TKey, TValue>
    {
        void Add(TKey key, IEnumerable<TValue> valueList);
        List<TValue> this[TKey key] { get; set;}
        bool Remove(TKey key, TValue value);
        void Add(TKey key, TValue value);
        bool ContainsKey(TKey key);

        ICollection<TKey> Keys {get;}
        bool Remove(TKey key);
        ICollection<TValue> Values{get;}

        void Add(KeyValuePair<TKey, TValue> item);
        void Clear();
        bool Contains(TKey key, TValue item);
        void CopyTo(KeyValuePair<TKey, TValue>[] array, int arrayIndex);
        int Count {get;}
        bool Remove(KeyValuePair<TKey, TValue> item);
    }
}
