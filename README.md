s2-geometry-library-csharp
==========================

This is a port of Google's S2 Geometry Library from both Java and C++
https://code.google.com/p/s2-geometry-library-java/

https://code.google.com/p/s2-geometry-library/

This library is can be used to create GeoHashes for fast querying. The Java version is used by AWS for 
GeoSpatial queries in DynamoDB.

S2 uses Hilbert Curves extensivly. 
For more info, see the original google presentation https://docs.google.com/presentation/d/1Hl4KapfAENAOf4gv-pSngKwvS_jwNVHRPZTTDzXXn6Q/view


Current status
---
Ready on NuGet, have fun!
`Install-Package S2Geometry`