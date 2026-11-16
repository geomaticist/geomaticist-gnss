# Sample data: NOAA CORS, GODE and GODN

One hour of real observations from two NOAA CORS stations at NASA Goddard (Greenbelt, Maryland, USA),
about 65 m apart, with the broadcast orbits of the day.

| file | what |
|---|---|
| `gode1660.25o.gz` | GODE (Septentrio PolaRx5TR, AOAD/M_T JPLA), 2025-06-15 12:00:00 to 12:59:30 GPS time, 30 s, RINEX 2.11 |
| `godn1660.25o.gz` | GODN (Javad TRE_3 Delta, TPSCR.G3 SCIS), same hour |
| `brdc1660.25n.gz` | GPS broadcast navigation, whole day |
| `brdc1660.25g.gz` | GLONASS broadcast navigation, whole day |

Try it:

```
gnss rinex info data/cors/gode1660.25o.gz
gnss baseline data/cors/gode1660.25o.gz data/cors/godn1660.25o.gz \
     --base-xyz 1130773.4883,-4831253.5703,3994200.4626 \
     --nav data/cors/brdc1660.25n.gz --nav data/cors/brdc1660.25g.gz
```

The base position is the approximate position from the GODE file header, which is fine for a demonstration.
Precise work also needs an ANTEX file (`--antex`, see the main README) and a published station position.

## Source and licence

Downloaded from the NOAA CORS Network archive (<https://geodesy.noaa.gov/CORS/>, files served from
<https://noaa-cors-pds.s3.amazonaws.com/rinex/2025/166/>). The daily station files were cut to one hour with
`gnss rinex trim` (which also leaves out the BeiDou and QZSS signals RINEX 2.11 cannot hold) and compressed
with gzip; the navigation files are as downloaded.

Data courtesy of the NOAA National Geodetic Survey CORS Network. As a work of the US Government it is in the
public domain in the United States.
