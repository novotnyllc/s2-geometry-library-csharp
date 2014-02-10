using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry.MultiMap;


namespace System.Collections.Generic
{
    internal class MultiMap<TKey, TValue> : IMultiMap<TKey,TValue>, IDictionary<TKey,TValue>
    {
        private Dictionary<TKey, List<TValue>> _interalStorage = new Dictionary<TKey, List<TValue>>();

        public MultiMap(){}

        public MultiMap(IEnumerable<KeyValuePair<TKey, TValue>> initialData)
        {
            foreach (var item in initialData)
            {
                Add(item);
            }
        }

        public void Add(TKey key, TValue value)
        {
            if (!_interalStorage.ContainsKey(key))
            {
                _interalStorage.Add(key, new List<TValue>());
            }
            _interalStorage[key].Add(value);
        }

        public void Add(TKey key, IEnumerable<TValue> valueList)
        {
            if (!_interalStorage.ContainsKey(key))
            {
                _interalStorage.Add(key, new List<TValue>());
            }
            foreach (TValue value in valueList)
            {
                _interalStorage[key].Add(value);
            }
        }

        public bool ContainsKey(TKey key)
        {
            return _interalStorage.ContainsKey(key);
        }

        public ICollection<TKey> Keys
        {
            get { return _interalStorage.Keys; }
        }

        public bool Remove(TKey key)
        {
            return _interalStorage.Remove(key);
        }

        bool IDictionary<TKey, TValue>.TryGetValue(TKey key, out TValue value)
        {
            if (!_interalStorage.ContainsKey(key))
            {
                value = default(TValue);
                return false;
            }
            value = _interalStorage[key].Last();
            return true;
        }

        public ICollection<TValue> Values
        {
            get 
            { 
                List<TValue> retVal = new List<TValue>();
                foreach (var item in _interalStorage)
                {
                    retVal.AddRange(item.Value);
                }
                return retVal;
            }
        }

        TValue IDictionary<TKey, TValue>.this[TKey key]
        {
            get
            {
                return _interalStorage[key].LastOrDefault();
            }
            set
            {
                Add(key,value);
            }
        }

        public void Add(KeyValuePair<TKey, TValue> item)
        {
            if (!_interalStorage.ContainsKey(item.Key))
            {
                _interalStorage.Add(item.Key, new List<TValue>());
            }
            _interalStorage[item.Key].Add(item.Value);
        }

        public void Clear()
        {
            _interalStorage.Clear();
        }

        public bool Contains(KeyValuePair<TKey, TValue> item)
        {
            List<TValue> valueList;
            if (_interalStorage.TryGetValue(item.Key, out valueList)) 
                return valueList.Contains(item.Value);
            return false;
        }

        public void CopyTo(KeyValuePair<TKey, TValue>[] array, int arrayIndex)
        {
            int i = arrayIndex;
            foreach (var item in _interalStorage)
            {
                foreach (TValue value in item.Value)
                {
                    array[i] = new KeyValuePair<TKey, TValue>(item.Key, value);
                    ++i;
                }
            }
        }

        public int Count
        {
            get
            {
                int count = 0;
                foreach (var item in _interalStorage)
                {
                    count += item.Value.Count;
                }
                return count;
            }
        }

        public bool CountIsAtLeast(int value)
        {
            int count = 0;
            foreach (var item in _interalStorage)
            {
                count += item.Value.Count;
                if (count >= value)
                    return true;
            }
            return false;
        }

        int ICollection<KeyValuePair<TKey,TValue>>.Count
        {
	        get { return _interalStorage.Count; }
        }

        bool ICollection<KeyValuePair<TKey, TValue>>.IsReadOnly
        {
            get { return false; }
        }

        public bool Remove(KeyValuePair<TKey, TValue> item)
        {
            if (!ContainsKey(item.Key)) return false;

            var list = _interalStorage[item.Key];
            var removed = list.Remove(item.Value);
            if (list.Count == 0)
                _interalStorage.Remove(item.Key); // clear out the dict
            
            return removed;
        }

        public IEnumerator<KeyValuePair<TKey, TValue>> GetEnumerator()
        {
            return new MultiMapEnumerator<TKey,TValue>(this);
        }

        public IEnumerable<KeyValuePair<TKey, TValue>>  SortedValues
        {
            get { return new SortedMultiMapEnumerable<TKey, TValue>(this); }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return new MultiMapEnumerator<TKey, TValue>(this);
        }


        public List<TValue> this[TKey key]
        {
            get
            {
                if (!_interalStorage.ContainsKey(key))
                    return new List<TValue>();
                return _interalStorage[key];
            }
            set
            {
                if (!_interalStorage.ContainsKey(key)) 
                    _interalStorage.Add(key, value);
                else _interalStorage[key] = value;
            }
        }

        public bool Remove(TKey key, TValue value)
        {
            if (!ContainsKey(key)) return false;
            return _interalStorage[key].Remove(value);
        }



        public bool Contains(TKey key, TValue item)
        {
           if (!_interalStorage.ContainsKey(key)) return false;
           return _interalStorage[key].Contains(item);
        }
    }
}
