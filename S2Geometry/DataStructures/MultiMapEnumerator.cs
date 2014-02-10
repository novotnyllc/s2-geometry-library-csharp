using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace System.Collections.Generic
{
    internal class MultiMapEnumerator<TKey,TValue> : IEnumerator<KeyValuePair<TKey,TValue>>
    {
        MultiMap<TKey,TValue> _map;
        IEnumerator<TKey> _keyEnumerator;
        IEnumerator<TValue> _valueEnumerator;

        public MultiMapEnumerator(MultiMap<TKey,TValue> map)
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
            _keyEnumerator = _map.Keys.GetEnumerator();
            _valueEnumerator = new List<TValue>().GetEnumerator();
        }
    }
}
