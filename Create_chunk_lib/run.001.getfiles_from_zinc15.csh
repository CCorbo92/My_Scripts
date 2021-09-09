#!/bin/csh -f 
### This script downloads 3D drug-like molecules from ZINC 15 (see readme_2019)


# Downloading mol2.gz files

set count = 0
foreach i (`cat ZINC-downloader-3D-mol2.gz.uri`)
    wget ${i}
    count ++
end
echo "Total number of tranches is: "
echo ${count}
exit
