## EC4MET 2.2

### (2025.10.07)

-   `plus.yr` parameter added to `get.SILO.weather()` function to download two years of weather data at a time for environments where estimated crop growth stages over-run the end of the sowing year.
-   `crop.locs()` added to easily define a grid of locations within the Australian crop growing regions.

## EC4MET 2.1

## EC4MET 2.0

### (2025.08.20)

-   Improved methods to bulk download soil tiff files to a local directory from the SLGA using the `dl.slga()` function. You can then point `get.S.ECs()` to this directory for faster processing.

## EC4MET 1.0

-   Initial development version.

-   SMI estimate function (`add.SMI()`) added.
