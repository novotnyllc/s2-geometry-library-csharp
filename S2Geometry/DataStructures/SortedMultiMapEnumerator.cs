using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry.MultiMap
{
    internal class SortedMultiMapEnumerable<TKey, TValue> : IEnumerable<KeyValuePair<TKey, TValue>>, IEnumerator<KeyValuePair<TKey, TValue>>
    {
        MultiMap<TKey, TValue> _map;
        IEnumerator<TKey> _keyEnumerator;
        IEnumerator<TValue> _valueEnumerator;


        public SortedMultiMapEnumerable(MultiMap<TKey, TValue> map)
        {
            this._map = map;
            Reset();
        }

        object IEnumerator.Current
        {
            get
            {
                return Current;
            }
        }

        public KeyValuePair<TKey,TValue> Current
        {
            get
            {
                return new KeyValuePair<TKey, TValue>(_keyEnumerator.Current, _valueEnumerator.Current);
            }
        }


        public void Dispose()
        {
            _keyEnumerator = null;
            _valueEnumerator = null;
            _map = null;
        }


        public bool MoveNext()
        {
            if (!_valueEnumerator.MoveNext())
            {
                if (!_keyEnumerator.MoveNext())
                    return false;
                _valueEnumerator = _map[_keyEnumerator.Current].GetEnumerator();
                _valueEnumerator.MoveNext();
                return true;
            }
            return true;
        }

        public void Reset()
        {
            _keyEnumerator = _map.Keys.OrderBy(k => k).GetEnumerator();
            _valueEnumerator = new List<TValue>().GetEnumerator();
        }

        public IEnumerator<KeyValuePair<TKey, TValue>> GetEnumerator()
        {
            return new SortedMultiMapEnumerable<TKey, TValue>(_map);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
