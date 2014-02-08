s2-geometry-library-csharp
==========================

This is a port of Google's S2 Geometry Library from both Java and C++
https://code.google.com/p/s2-geometry-library-java/
https://code.google.com/p/s2-geometry-library/


Current status
---
Two failing tests:
- testNeighbors (on level 29)
- testRandomCaps (possibly related the the other failing test)

Code was ported from Java/C++. It still needs to be cleaned up to make it C# friendly.
The current effort is ensuring proper functionality (tests all pass).