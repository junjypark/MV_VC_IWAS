
## Data obtaining
**1. [GWAS summary statistics on 3,935 Imaging Derived Phenotypes (IDPs)](https://open.win.ox.ac.uk/ukbiobank/big40/)**

1.1 **download** located in: /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP
  for i in {0001..3935}; do
    curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/${i}.txt.gz
  done
```

1.2 **unzip**
```
  for i in {0001..3935}; do
    gunzip ${i}.txt.gz
  done
```

1.3 **separate to chromosome**, e.g. chr19
```
  for i in {0001..3935}; do
    awk -F' ' 'NR==1{print;next}$1 == 19' /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP/${i}.txt > /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/IDP${i}_chr19.txt
  done
```

**2. [phase3 1000G reference panel, build 37, 2504 samples](https://www.cog-genomics.org/plink/2.0/resources)**

2.1 **download 1000G** 
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/rawdata
  wget -O all_phase3.psam "https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1"
  wget -O all_phase3.pgen.zst "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
  wget -O all_phase3.pvar.zst "https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"
```

2.2 **unzip**
```
  /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
  /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
```

**3. [Stage 1 International Genomics of Alzheimer’s Project (IGAP) ](https://www.niagads.org/datasets/ng00036)**

3.1 **download IGAP** 
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/IGAP
  wget -O IGAP_summary_statistics.zip     
  "https://www.niagads.org/system/tdf/public_docs/IGAP_summary_statistics.zip?file=1&type=field_collect on_item&id=185&force="
```

3.2 **unzip**
```
  gunzip IGAP_summary_statistics.zip
```

## Data filterting

### UKB
Filter for EUR in later steps

### 1000G

**4. Choose EUR population, 263 females, 240 males**

4.1 **Prepare sub-population filter file**
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/rawdata
  awk 'NR == 1 || $5 == "EUR" {print $1}' all_phase3.psam > EUR_1kg_samples.txt
```

4.2 **Generate sub-population fileset**
```
  /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --pfile all_phase3 --keep EUR_1kg_samples.txt --make-pgen --out EUR_phase3_autosomes
```

**5. Remove duplicated variants**
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/rawdata
  /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --pfile EUR_phase3_autosomes --rm-dup force-first --make-pgen --out EUR_phase3_DuplicatesRemoved
```

**6. Separate 1000G data to chr**

e.g. chr19
```
  for i in {1..22}; do
    /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --pfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/rawdata/EUR_phase3_DuplicatesRemoved --allow-extra-chr --chr ${i} --make-pgen --out /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR/1000G_all_phase3_EUR_DuplicatesRemoved_chr${i}
  done
```

**7. LD clumping with window 500K and $r^2$ 0.5**

e.g. chr22
```
  idp_ids="1452 1453 1454 1455 1456 1457 1458 1459 1460 1461 1462 1463 1464 1465 1466 1467 1468 1469 1470 1471 1472 1473 1474 1475 1476 1477 1478 1479 1480 1481 1482 1483 1484 1485 1486 1487 1488 1489 1490 1491 1492 1493 1494 1495 1496 1497 1498 1499 1500 1501 1502 1503 1504 1505 1506 1507 1508 1509 1510 1511 1512 1513 1514 1515 1516 1517 1518 1519 1520 1521 1522 1523 1524 1525 1526 1527 1528 1529 1530 1531 1532 1533 1534 1535 1536 1537 1538 1539 1540 1541 1542 1543 1544 1545 1546 1547 1548 1549 1550 1551 1552 1553 1554 1555 1556 1557 1558 1559 1560 1561 1562 1563 1564 1565 1566 1567 1568 1569 1570 1571 1572 1573 1574 1575 1576 1577 1578 1579 1580 1581 1582 1583 1584 1585 1586 1587 1588 1589 1590 1591 1592 1593 1594 1595 1596 1597 1598 1599 1600 1601 1602 1603 1604 1605 1606 1607 1608 1609 1610 1611 1612 1613 1614 1615 1616 1617 1618 1619 1620 1621 1622 1623 1624 1625 1626 1627 1628 1629 1630 1631 1632 1633 1634 1635 1636 1637 1638 1639 1640 1641 1642 1643 1644 1645 1646 1647 1648 1649 1650 1651 1652 1653 1654 1655 1656 1657 1658 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671 1672 1673 1674 1675 1676 1677 1678 1679 1680 1681 1682 1683 1684 1685 1686 1687 1688 1689 1690 1691 1692 1693 1694 1695 1696 1697 1698 1699 1700 1701 1702 1703 1704 1705 1706 1707 1708 1709 1710 1711 1712 1713 1714 1715 1716 1717 1718 1719 1720 1721 1722 1723 1724 1725 1726 1727 1728 1729 1730 1731 1732 1733 1734 1735 1736 1737 1738 1739 1740 1741 1742 1743 1744 1745 1746 1747 1748 1749 1750 1751 1752 1753 1754 1755 1756 1757 1758 1759 1760 1761 1762 1763 1764 1765 1766 1767 1768 1769 1770 1771 1772 1773 1774 1775 1776 1777 1778 1779 1780 1781 1782 1783 1784 1785 1786 1787 1788 1789 1790 1791 1792 1793 1794 1795 1796 1797 1798 1799 1800 1801 1802 1803 1804 1805 1806 1807 1808 1809 1810 1811 1812 1813 1814 1815 1816 1817 1818 1819 1820 1821 1822 1823 1824 1825 1826 1827 1828 1829 1830 1831 1832 1833 1834 1835 1836 1837 1838 1839 1840 1841 1842 1843 1844 1845 1846 1847 1848 1849 1850 1851 1852 1853 1854 1855 1856 1857 1858 1859 1860 1861 1862 1863 1864 1865 1866 1867 1868 1869 1870 1871 1872 1873 1874 1875 1876 1877 1878 1879 1880 1881 1882 1883 1884 1885 1886 1887 1888 1889 1890 1891 1892 1893 1894 1895 1896 1897 1898 1899 1900 1901 1902 1903 1904 1905 1906 1907 1908 1909 1910 1911 1912 1913 1914 1915 1916 1917 1918 1919 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040 2041 2042 2043 2044 2045 2046 2047 2048 2049 2050 2051 2052 2053 2054 2055 2056 2057 2058 2059 2060 2061 2062 2063 2064 2065 2066 2067 2068 2069 2070 2071 2072 2073 2074 2075 2076 2077 2078 2079 2080 2081 2082 2083 2084 2085 2086 2087 2088 2089 2090 2091 2092 2093 2094 2095 2096 2097 2098 2099 2100 2101 2102 2103 2104 2105 2106 2107 2108 2109 2110 2111 2112 2113 2114 2115 2116 2117 2118 2119 2120 2121 2122 2123 2124 2125 2126 0011 0012 0013 0014 0015 0016 0017 0018 0019 0020 0021 0022 0023 0024 0026 0027 0028 0029 0030 0031 0032 0033 0034 0035 0036 0037 0038 0039 0040 0041 0042 0043 0044 0045 0046 0047 0048 0049 0050 0051 0052 0053 0054 0055 0056 0057 0058 0059 0060 0061 0062 0063 0064 0065 0066 0067 0068 0069 0070 0071 0072 0073 0074 0075 0076 0077 0078 0079 0080 0081 0082 0083 0084 0085 0086 0087 0088 0089 0090 0091 0092 0093 0094 0095 0096 0097 0098 0099 0100 0101 0102 0103 0104 0105 0106 0107 0108 0109 0110 0111 0112 0113 0114 0115 0116 0117 0118 0119 0120 0121 0122 0123 0124 0125 0126 0127 0128 0129 0130 0131 0132 0133 0134 0135 0136 0137 0138 0139 0140 0141 0142 0143 0144 0145 0146 0147 0148 0149 0150 0151 0152 0153 0154 0155 0156 0157 0158 0159 0160 0161 0162 0163 0164 0165 0166 0167 0168 0169 0170 0171 0172 0173 0174 0175 0176 0177 0178 0179 0180 0181 0182 0183 0184 0185 0186 0187 0188 0189 0190 0191 0192 0193 0194 0195 0196 0197 0198 0199 0200 0201 0202 0203 0204 0205 0206 0207 0208 0209 0210 0211 0212 0213 0214 0215 0216 0217 0218 0219 0220 0221 0222 0223 0224 0225 0226 0227 0228 0229 0230 0231 0232 0233 0234 0235 0236 0237 0238 0239 0240 0241 0242 0243 0244 0245 0246 0247 0248 0249 0250 0251 0252 0253 0254 0255 0256 0257 0258 0259 0260 0261 0262 0263 0264 0265 0266 0267 0268 0269 0270 0271 0272 0273 0274 0275 0276 0277 0278 0279 0280 0281 0282 0283 0284 0285 0286 0287 0288 0289 0290 0291 0292 0293 0294 0295 0296 0297 0298 0299 0300 0301 0302 0303 0304 0305 0306 0307 0308 0309 0310 0311 0312 0313 0314 0315 0316 0317 0318 0319 0320 0321 0322 0323 0324 0325 0326 0327 0328 0329 0330 0331 0332 0333 0334 0335 0336 0337 0338 0339 0340 0341 0342 0343 0344 0345 0346 0347 0348 0349 0350 0351 0352 0353 0354 0355 0356 0357 0358 0359 0360 0361 0362 0363 0364 0365 0366 0367 0368 0369 0370 0371 0372 0373 0374 0375 0376 0377 0378 0379 0380 0381 0382 0383 0384 0385 0386 0387 0388 0389 0390 0391 0392 0393 0394 0395 0396 0397 0398 0399 0400 0401 0402 0403 0404 0405 0406 0407 0408 0409 0410 0411 0412 0413 0414 0415 0416 0417 0418 0419 0420 0421 0422 0423 0424 0425 0426 0427 0428 0429 0430 0431 0432 0433 0434 0435 0436 0437 0438 0439 0440 0441 0442 0443 0444 0445 0446 0447 0448 0449 0450 0451 0452 0453 0454 0455 0456 0457 0458 0459 0460 0461 0462 0463 0464 0465 0466 0467 0468 0469 0470 0471 0472 0473 0474 0475 0476 0477 0478 0479 0480 0481 0482 0483 0484 0485 0486 0487 0488 0489 0490 0491 0492 0493 0494 0495 0496 0497 0498 0499 0500 0501 0502 0503 0504 0505 0506 0507 0508 0509 0510 0511 0512 0513 0514 0515 0516 0517 0518 0519 0520 0521 0522 0523 0524 0525 0526 0527 0528 0529 0530 0531 0532 0533 0534 0535 0536 0537 0538 0539 0540 0541 0542 0543 0544 0545 0546 0547 0548 0549 0550 0551 0552 0553 0554 0555 0556 0557 0558 0559 0560 0561 0562 0563 0564 0565 0566 0567 0568 0569 0570 0571 0572 0573 0574 0575 0576 0577 0578 0579 0580 0581 0582 0583 0584 0585 0586 0587 0588 0589 0590 0591 0592 0593 0594 0595 0596 0597 0598 0599 0600 0601 0602 0603 0604 0605 0606 0607 0608 0609 0610 0611 0612 0613 0614 0615 0616 0617 0618 0619 0620 0621 0622 0623 0624 0625 0626 0627 0628 0629 0630 0631 0632 0633 0634 0635 0636 0637 0638 0639 0640 0641 0642 0643 0644 0645 0646 0647 0648 0649 0650 0651 0652 0653 0654 0655 0656 0657 0658 0659 0660 0661 0662 0663 0664 0665 0666 0667 0668 0669 0670 0671 0672 0673 0674 0675 0676 0677 0678 0679 0680 0681 0682 0683 0684 0685 0686 0687 0688 0689 0690 0691 0692 0693 0694 0695 0696 0697 0698 0699 0700 0701 0702 0703 0704 0705 0706 0707 0708 0709 0710 0711 0712 0713 0714 0715 0716 0717 0718 0719 0720 0721 0722 0723 0724 0725 0726 0727 0728 0729 0730 0731 0732 0733 0734 0735 0736 0737 0738 0739 0740 0741 0742 0743 0744 0745 0746 0747 0748 0749 0750 0751 0752 0753 0754 0755 0756 0757 0758 0759 0760 0761 0762 0763 0764 0765 0766 0767 0768 0769 0770 0771 0772 0773 0774 0775 0776 0777 0778 0779 0780 0781 0782 0783 0784 0785 0786 0787 0788 0789 0790 0791 0792 0793 0794 0795 0796 0797 0798 0799 0800 0801 0802 0803 0804 0805 0806 0807 0808 0809 0810 0811 0812 0813 0814 0815 0816 0817 0818 0819 0820 0821 0822 0823 0824 0825 0826 0827 0828 0829 0830 0831 0832 0833 0834 0835 0836 0837 0838 0839 0840 0841 0842 0843 0844 0845 0846 0847 0848 0849 0850 0851 0852 0853 0854 0855 0856 0857 0858 0859 0860 0861 0862 0863 0864 0865 0866 0867 0868 0869 0870 0871 0872 0873 0874 0875 0876 0877 0878 0879 0880 0881 0882 0883 0884 0885 0886 0887 0888 0889 0890 0891 0892 0893 0894 0895 0896 0897 0898 0899 0900 0901 0902 0903 0904 0905 0906 0907 0908 0909 0910 0911 0912 0913 0914 0915 0916 0917 0918 0919 0920 0921 0922 0923 0924 0925 0926 0927 0928 0929 0930 0931 0932 0933 0934 0935 0936 0937 0938 0939 0940 0941 0942 0943 0944 0945 0946 0947 0948 0949 0950 0951 0952 0953 0954 0955 0956 0957 0958 0959 0960 0961 0962 0963 0964 0965 0966 0967 0968 0969 0970 0971 0972 0973 0974 0975 0976 0977 0978 0979 0980 0981 0982 0983 0984 0985 0986 0987 0988 0989 0990 0991 0992 0993 0994 0995 0996 0997 0998 0999 1000 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020 1021 1022 1023 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 1035 1036 1037 1038 1039 1040 1041 1042 1043 1044 1045 1046 1047 1048 1049 1050 1051 1052 1053 1054 1055 1056 1057 1058 1059 1060 1061 1062 1063 1064 1065 1066 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 1077 1078 1079 1080 1081 1082 1083 1084 1085 1086 1087 1088 1089 1090 1091 1092 1093 1094 1095 1096 1097 1098 1099 1100 1101 1102 1103 1104 1105 1106 1107 1108 1109 1110 1111 1112 1113 1114 1115 1116 1117 1118 1119 1120 1121 1122 1123 1124 1125 1126 1127 1128 1129 1130 1131 1132 1133 1134 1135 1136 1137 1138 1139 1140 1141 1142 1143 1144 1145 1146 1147 1148 1149 1150 1151 1152 1153 1154 1155 1156 1157 1158 1159 1160 1161 1162 1163 1164 1165 1166 1167 1168 1169 1170 1171 1172 1173 1174 1175 1176 1177 1178 1179 1180 1181 1182 1183 1184 1185 1186 1187 1188 1189 1190 1191 1192 1193 1194 1195 1196 1197 1198 1199 1200 1201 1202 1203 1204 1205 1206 1207 1208 1209 1210 1211 1212 1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 1225 1226 1227 1228 1229 1230 1231 1232 1233 1234 1235 1236 1237 1238 1239 1240 1241 1242 1243 1244 1245 1246 1247 1248 1249 1250 1251 1252 1253 1254 1255 1256 1257 1258 1259 1260 1261 1262 1263 1264 1265 1266 1267 1268 1269 1270 1271 1272 1273 1274 1275 1276 1277 1278 1279 1280 1281 1282 1283 1284 1285 1286 1287 1288 1289 1290 1291 1292 1293 1294 1295 1296 1297 1298 1299 1300 1301 1302 1303 1304 1305 1306 1307 1308 1309 1310 1311 1312 1313 1314 1315 1316 1317 1318 1319 1320 1321 1322 1323 1324 1325 1326 1327 1328 1329 1330 1331 1332 1333 1334 1335 1336 1337 1338 1339 1340 1341 1342 1343 1344 1345 1346 1347 1348 1349 1350 1351 1352 1353 1354 1355 1356 1357 1358 1359 1360 1361 1362 1363 1364 1365 1366 1367 1368 1369 1370 1371 1372 1373 1374 1375 1376 1377 1378 1379 1380 1381 1382 1383 1384 1385 1386 1387 1388 1389 1390 1391 1392 1393 1394 1395 1396 1397 1398 1399 1400 1401 1402 1403 1404 1405 1406 1407 1408 1409 1410 1411 1412 1413 1414 1415 1416 1417 1418 1419 1420 1421 1422 1423 1424 1425 1426 1427 1428 1429 1430 1431 1432 1433 1434 1435 1436 2219 2220 2221 2222 2223 2224 2225 2226 2227 2228 2229 2230 2231 2232 2233 2234 2235 2236 2237 2238 2239 2240 2241 2242 2243 2244 2245 2246 2247 2248 2249 2250 2251 2252 2253 2254 2255 2256 2257 2258 2259 2260 2261 2262 2263 2264 2265 2266 2267 2268 2269 2270 2271 2272 2273 2274 2275 2276 2277 2278 2279 2280 2281 2282 2283 2284 2285 2286 2287 2288 2289 2290 2291 2292 2293 2294 2295 2296 2297 2298 2299 2300 2301 2302 2303 2304 2305 2306 2307 2308 2309 2310 2311 2312 2313 2314 2315 2316 2317 2318 2319 2320 2321 2322 2323 2324 2325 2326 2327 2328 2329 2330 2331 2332 2333 2334 2335 2336 2337 2338 2339 2340 2341 2342 2343 2344 2345 2346 2347 2348 2349 2350 2351 2352 2353 2354 2355 2356 2357 2358 2359 2360 2361 2362 2363 2364 2365 2366 2367 2368 2369 2370 2371 2372 2373 2374 2375 2376 2377 2378 2379 2380 2381 2382 2383 2384 2385 2386 2387 2388 2389 2390 2391 2392 2393 2394 2395 2396 2397 2398 2399 2400 2401 2402 2403 2404 2405 2406 2407 2408 2409 2410 2411 2412 2413 2414 2415 2416 2417 2418 2419 2420 2421 2422 2423 2424 2425 2426 2427 2428"
```
```
  for id in $idp_ids; do
    idp_input="/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/IDP${id}_chr22.txt"
    idp_output="/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/clump_0.5/chr22/IDP${id}_chr22_0.5"
    /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --pfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR/1000G_all_phase3_EUR_DuplicatesRemoved_chr22 --clump "$idp_input" --clump-id-field rsid --clump-log10 --clump-p-field "pval(-log10)" --clump-a1-field a1 --clump-r2 0.50 --clump-kb 500 --out "$idp_output"
  done
```

**8. Obtain clumped IDP SNPs list**
```
  IDPSNPs_func <- function(chr) {
    library(data.table)
    all_MRI_ids <- readRDS("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/all_MRI_ids.rds")
    for (IDP in all_MRI_ids) {
      tryCatch({
        formatted_IDP <- sprintf("%04d", IDP)
        filename <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/clump_0.5/chr", chr, "/IDP", formatted_IDP, "_chr", chr, "_0.5.clumps", sep = "")
        clumped <- fread(filename)
        clumped_rs <- clumped$ID
        clumped_rs <- as.data.frame(clumped_rs)
        outfilename <- paste("/chr", chr, "/SNPlist_IDP",  formatted_IDP, "_chr", chr, "_0.5.txt", sep = "")
        write.table(clumped_rs, outfilename, row.names = F, col.names = F, quote = F)
      }, error = function(e){print(IDP)})
    }
  }
```
```
  for (i in 1:22) {
    IDPSNPs_func(chr = i)
  }
```

### IGAP

**9. IGAP split by chr**
```
  library(data.table)
  summary_data <- fread("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_stage_1.txt")
  for (i in 1:22) {
    dft <- summary_data[which(summary_data$Chromosome == i), ]
    filename <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_chr", i, ".txt", sep = "")
    write.table(dft, file = filename)
  }
```

**10. IGAP data split by gene**
```
  library(data.table)
  split_gene <- function(chr) {
    # IGAP summary data
    filename_IGAP <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_chr", chr, ".txt", sep = "")
    summary_data_chr <- fread(filename_IGAP)
  
    # get max gene numbers
    path_value = paste("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/genelist1M/chr", chr, sep = "")
    filenames_genenum <- list.files(path=path_value, full.names = F)
    gene_n <- max(as.integer(sub(".*gene(\\d+)\\.txt", "\\1", filenames_genenum)))
    print(gene_n)
  
    # get IGAP by gene
    for (i in 1:gene_n) {
      tryCatch({
        filename <- paste("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/genelist1M/chr", chr, "/hg19_chr", chr, "_gene", i, ".txt", sep = "")
        snprange <- read.table(file = filename)
        listofsnps <- summary_data_chr[which(summary_data_chr$Position >= snprange$V2 & summary_data_chr$Position <= snprange$V3), ]
        listofsnps <- as.data.frame(listofsnps)
        outfilename <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/genelist/chr", chr, "/IGAP_genelist_chr", chr, "_gene", i, ".txt", sep = "")
        write.table(listofsnps, file = outfilename, row.names = F, col.names = F, quote = F)
      }, error = function(e){print(i)})
    }
  }

  for (chr in 1:22) {
    split_gene(chr = chr)
    print(chr)
  }
```

### Intersect 1000G, UKB, and IGAP SNPs

**11. Intersect clumped 1000G, IGAP and UKB by chr, IDP, and gene**

Use *intersection_0.5.R* and *intersection_0.5.sh* in

*/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5*
```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5
  for set in {1..168}; do
    sbatch --export=set=$set,name="intersection_0.5" intersection_0.5.sh -o intersection_0.5_result_${set}.out
  done
```

**12. Summarize SNPs overlap between modalities by chr by gene**

12.1 **SNPs in each gene by IDP and modility**

```
  library(data.table)
  checkSNP_IDP_gene <- function(chr, geneindex) {
    tryCatch({
  
    MRI_IDS <- readRDS("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.2/MRI_IDS.rds")
  
    chr_path <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5/chr", chr, "/", sep = "")
    all_files <- list.files(chr_path)
    matching_files <- grep(paste0("gene", geneindex, ".txt"), all_files, value = TRUE)
  
    N_snp <- rep(NA, length(matching_files))
    SNPs <- rep(NA, length(matching_files))
  
    IDP_pattern <- paste0("Intersect3_0.5_chr", chr, "_IDP(\\d+)_gene\\d+\\.txt")
    IDP_ID <- sub(IDP_pattern, "\\1", matching_files)
  
    IDP_type <- rep(NA, length(matching_files))
  
    for (IDPindex in 1:length(matching_files)) {
      filename <- paste0(chr_path, matching_files[IDPindex])
      SNPlist <- fread(filename, header = F)
      N_snp[IDPindex] <- nrow(SNPlist)
      SNPs[IDPindex] <- paste(unlist(SNPlist[, 1]), collapse = ", ", sep = "")
      if (as.numeric(IDP_ID[IDPindex]) %in% MRI_IDS$dMRI_ids) {
        IDP_type[IDPindex] <- "D"
      } else if (as.numeric(IDP_ID[IDPindex]) %in% MRI_IDS$sMRI_ids) {
        IDP_type[IDPindex] <- "S"
      } else if (as.numeric(IDP_ID[IDPindex]) %in% MRI_IDS$fMRI_ids) {
        IDP_type[IDPindex] <- "F"
      }
    }
  
    SNPs_IDP_gene <- data.frame(IDP_ID, IDP_type, N_snp, SNPs)
  
    outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5_bygene/chr", chr, "/SNPresults_0.5_chr", chr, "_gene", geneindex, ".txt")
    write.table(SNPs_IDP_gene, outname, col.names = T, row.names = F)
  
    }, error = function(e){print(geneindex)})
  }
```

12.1 **Overlapping by chr**

```
  library(data.table)
    overlap <- function(chr) {
  
    gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
  
    chrresults <- data.frame(matrix(NA, nrow = gene_n_v[chr], ncol = 7))
    colnames(chrresults) <- c("geneindex", "nSNPs_S", "nSNPs_D", "nSNPs_F", "SandD", "SandF", "DandF")
    chrresults[, 1] <- 1:gene_n_v[chr]
  
    for (geneindex in 1:gene_n_v[chr]) {
      tryCatch({
      
        readname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5_bygene/chr", chr, 
                           "/SNPresults_0.5_chr", chr, "_gene", geneindex, ".txt")
        readdata <- data.frame(fread(readname, header = T))
      
        SNPs_S <- unlist(readdata$SNPs[which(readdata$IDP_type == "S")])
        SNPs_S <- unique(unlist(strsplit(SNPs_S, ", ")))
        chrresults[geneindex, "nSNPs_S"] <- length(SNPs_S)
        SNPs_D <- unlist(readdata$SNPs[which(readdata$IDP_type == "D")])
        SNPs_D <- unique(unlist(strsplit(SNPs_D, ", ")))
        chrresults[geneindex, "nSNPs_D"] <- length(SNPs_D)
        SNPs_F <- unlist(readdata$SNPs[which(readdata$IDP_type == "F")])
        SNPs_F <- unique(unlist(strsplit(SNPs_F, ", ")))
        chrresults[geneindex, "nSNPs_F"] <- length(SNPs_F)
      
        if ((length(SNPs_S) > 0) & (length(SNPs_D) > 0)) {
          chrresults[geneindex, "SandD"] <- length(which(SNPs_S %in% SNPs_D))
        }
        if ((length(SNPs_D) > 0) & (length(SNPs_F) > 0)) {
          chrresults[geneindex, "DandF"] <- length(which(SNPs_D %in% SNPs_F))
        }
        if ((length(SNPs_S) > 0) & (length(SNPs_F) > 0)) {
          chrresults[geneindex, "SandF"] <- length(which(SNPs_S %in% SNPs_F))
        }
      
      }, error = function(e){print(geneindex)})
    }
    outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_overlap_0.5_bychr/overlap_chr", chr, ".txt")
    write.table(chrresults, outname, row.names = F, col.names = T)
  }

  for (chr in 1:22) {
    overlap(chr)
  }
```

**13. Summarize how many genes in the chr shared at least 1 SNPs in at least 2 modalities**

```
  atleast2 <- data.frame(matrix(NA, nrow = 22, ncol = 3))
  colnames(atleast2) <- c("chr", "size", "overlap")
  atleast2[, 1] <- 1:22
  gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
  atleast2[, 2] <- gene_n_v
  for (chr in 1:22) {
    readname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_overlap_0.5_bychr/overlap_chr", chr, ".txt")
    chrresults <- data.frame(fread(readname, header = T))
    chroverlap <- rep(NA, gene_n_v[chr])
    for (geneindex in 1:gene_n_v[chr]) {
      chroverlap[geneindex] <- !all(unlist(chrresults[geneindex, 5:7]) %in% c(NA, 0))
    }
    atleast2[chr, 3] <- sum(chroverlap)
  }

  write.table(atleast2, "/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_overlap_0.5_bychr/summary.txt", row.names = F, col.names = T)
```

**14. SNPs list after intersection for LD**

 14.1 **By gene:**

```
  library(data.table)
  snplist <- function(chr) {
  
    gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
  
    for (geneindex in 1:gene_n_v[chr]) {
      tryCatch({
      
        readname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5_bygene/chr", chr, 
                           "/SNPresults_0.5_chr", chr, "_gene", geneindex, ".txt")
        readdata <- data.frame(fread(readname, header = T))
      
        SNPs_S <- unlist(readdata$SNPs[which(readdata$IDP_type == "S")])
        SNPs_S <- unique(unlist(strsplit(SNPs_S, ", ")))
        SNPs_D <- unlist(readdata$SNPs[which(readdata$IDP_type == "D")])
        SNPs_D <- unique(unlist(strsplit(SNPs_D, ", ")))
        SNPs_F <- unlist(readdata$SNPs[which(readdata$IDP_type == "F")])
        SNPs_F <- unique(unlist(strsplit(SNPs_F, ", ")))
      
        SNPs <- unique(c(SNPs_S, SNPs_D, SNPs_F))
        SNPs <- data.frame(SNPs)
      
        outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5/chr", chr, 
                        "/SNPslist_0.5_chr", chr, "_gene", geneindex, ".txt")
        write.table(SNPs, outname, row.names = F, col.names = F, quote = FALSE)
      
      }, error = function(e){print(geneindex)})
    }
  }

  for (chr in 1:22) {
    snplist(chr)
  }
```

14.2 **By chr:**

```
  library(data.table)
  snplist_chr <- function(chr) {
  
    gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
    SNP <- c()
  
    for (geneindex in 1:gene_n_v[chr]) {
      tryCatch({
      
        readname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5_bygene/chr", chr, 
                           "/SNPresults_0.5_chr", chr, "_gene", geneindex, ".txt")
        readdata <- data.frame(fread(readname, header = T))
      
        SNPs_S <- unlist(readdata$SNPs[which(readdata$IDP_type == "S")])
        SNPs_S <- unique(unlist(strsplit(SNPs_S, ", ")))
        SNPs_D <- unlist(readdata$SNPs[which(readdata$IDP_type == "D")])
        SNPs_D <- unique(unlist(strsplit(SNPs_D, ", ")))
        SNPs_F <- unlist(readdata$SNPs[which(readdata$IDP_type == "F")])
        SNPs_F <- unique(unlist(strsplit(SNPs_F, ", ")))
      
        SNPs <- unique(c(SNPs_S, SNPs_D, SNPs_F))
        SNP <- unlist(c(SNP, SNPs))
      
      }, error = function(e){print(geneindex)})
    
      SNP <- unique(SNP)
      outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5", 
                        "/SNPslist_0.5_chr", chr, ".txt")
      write.table(SNP, outname, row.names = F, col.names = F, quote = FALSE)
    
    }
  }

  for (chr in 1:22) {
    snplist_chr(chr)
  }
```

## Prepare data

### Prepare UKB summary data

**15. Calculate UKB MAF by chr**
```
  for chr in {1..22}; do
    /gpfs/fs1/home/j/junpark/tianyu47/plink2/plink2 --bfile /gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/UKBB_ImgSampleGenotype_Chr${chr} --freq --out /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKB_MAF_bychr/UKBB_ImgSampleGenotype_MAF_Chr${chr}
  done
```

**16. UKB summary data by chr, gene, IDP**

Use *UKBGWAS_0.5.R* and *UKBGWAS_0.5.sh* in

*/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5*

```
  cd /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5
  for set in {1..168}; do
    sbatch --export=set=$set,name="UKBGWAS_0.5" UKBGWAS_0.5.sh -o UKBGWAS_0.5_result_${set}.out
  done
```

### 1000G LD matrix

**17. Prepare 1000G bfile by gene for LD calculation**

17.1 **filter 1000G for intersection snps by chr**

```
  for chr in {1..22}; do
    /gpfs/fs1/home/j/junpark/tianyu47/plink2_Nov/plink2 --pfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR/1000G_all_phase3_EUR_DuplicatesRemoved_chr${chr} --extract /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5/SNPslist_0.5_chr${chr}.txt --make-pgen --out /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR_filter/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr}
  done
```

17.2 **Convert pfile to bfile**

```
  for chr in {1..22}; do
    /gpfs/fs1/home/j/junpark/tianyu47/plink2_Nov/plink2 --pfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR_filter/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr} --make-bed --out /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR_bfile/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr}
  done
```

17.3 **Separate bfile data by gene**

```
  GENES_PER_CHR=(2568 1648 1359 1018 1172 1282 1235 919 1028 990 1523 1282 550 843 954 1083 1489 383 1806 723 330 583)
  for chr in {1..22}; do
    for ((geneindex=1; geneindex <= GENES_PER_CHR[chr-1]; geneindex++)); do
      plink --bfile /gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/UKBB_ImgSampleGenotype_Chr${chr} --extract /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5/chr${chr}/SNPslist_0.5_chr${chr}_gene${geneindex}.txt --make-bed --out /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBind_0.5/chr${chr}/UKBind_chr${chr}_gene${geneindex}
    done
  done
```

**18. LD calculation**

18.1 **LD using plink**
```
  GENES_PER_CHR=(2568 1648 1359 1018 1172 1282 1235 919 1028 990 1523 1282 550 843 954 1083 1489 383 1806 723 330 583)
  for chr in {1..22}; do
    for ((geneindex=1; geneindex <= GENES_PER_CHR[chr-1]; geneindex++)); do
      plink --bfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR_bfile/chr${chr}/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr}_gene${geneindex} --r square --write-snplist --out /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/LD_0.5/chr${chr}/LD_chr${chr}_gene${geneindex}
    done
  done
```

18.2 **Adding snp name**

```
  library(data.table)
  LD_snpname <- function(chr) {
    gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
    for (geneindex in 1:gene_n_v[chr]) {
      tryCatch({
        ldfilename <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".ld")
        ld <- as.data.frame(fread(ldfilename, header = F))
        snpnamelist <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".snplist")
        snplist <- unlist(fread(snpnamelist, header = F))
        colnames(ld) <- snplist
        rownames(ld) <- snplist
        outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/R_LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".rds")
        saveRDS(ld, outname)
      }, error = function(e){print(geneindex)})
    }
  }

  for (chr in 1:22) {
    LD_snpname(chr)
  }
```

### IGAP data

**19. IGAP data by chr by gene**

```
  library(data.table)
  IGAP_data <- function(chr) {
  
    gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
    IGAPfilename <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_chr", chr, ".txt")
    IGAP <- as.data.frame(fread(IGAPfilename))[, -1]
  
    for (geneindex in 1:gene_n_v[chr]) {
    
      tryCatch({
        snpnamelist <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".snplist")
        snplist <- unlist(fread(snpnamelist, header = F))
        IGAP_filtered <- IGAP[IGAP$MarkerName %in% snplist, ]
        IGAP_sorted <- IGAP_filtered[match(snplist, IGAP_filtered$MarkerName), ]
        outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_0.5_bygene/chr", chr, "/IGAP_chr", chr, "_gene", geneindex, ".rds")
        saveRDS(IGAP_sorted, outname)
      }, error = function(e){print(geneindex)})
    }
  }

  for (chr in 1:22) {
    IGAP_data(chr)
  }
```

### UKB individual data

**20. Filter UKB individual data in case needed**

```
  GENES_PER_CHR=(2568 1648 1359 1018 1172 1282 1235 919 1028 990 1523 1282 550 843 954 1083 1489 383 1806 723 330 583)
  for chr in {1..22}; do
    for ((geneindex=1; geneindex <= GENES_PER_CHR[chr-1]; geneindex++)); do
      plink --bfile /gpfs/fs0/scratch/j/junpark/tianyu47/1000G/chrdata_EUR_bfile/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr} --extract /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5/chr${chr}/SNPslist_0.5_chr${chr}_gene${geneindex}.txt --make-bed --out /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBind_0.5/chr${chr}/1000G_all_phase3_EUR_DuplicatesRemoved_filtered_chr${chr}_gene${geneindex}
    done
  done
```

### UKB GWAS

**21. Calculate GWAS statistics for all SNPs in UKB by gene**
Use *UKBGWAS.R*

```
library(dplyr)
library(plink2R)

GWASfunc <- function(chr, index) {
  tryCatch({
    a_name <- paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/genedata1M/chr", chr, "/chr", 
                     chr, "_gene", index)
    a <- read_plink(a_name)
    b <- readRDS("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/behavioral.rds")
    sort_b <- b %>%
      filter(eid %in% a$fam[, 2]) %>%
      arrange(match(eid, a$fam[, 2]))
    
    GWAS_result <- as.data.frame(matrix(NA, ncol(a$bed), 5))
    colnames(GWAS_result) <- c("rsID", "beta", "SE", "z", "p")
    
    for (snp in 1:ncol(a$bed)) {
      dat <- sort_b[, 1:14]
      dat <- cbind(dat, a$bed[, snp])
      names(dat)[15] <- colnames(a$bed)[snp]
      dat <- dat[, -2]
      dat$sex <- as.factor(dat$sex)
      fit <- glm(AD ~ ., data = dat, family = binomial(link = "logit"))
      GWAS_result[snp, 1] <- colnames(a$bed)[snp]
      GWAS_result[snp, 2:5] <- summary(fit)$coef[14, ]
    }
    
    outnames <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_rawADGWAS/chr", chr, "/UKB_GWAS_chr", chr, "_gene", index, ".rds")
    saveRDS(GWAS_result, outnames)
    
  }, error = function(e){print(index)})
}

gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)

for (chr in 1:22) {
  for (index in 1:gene_n_v[chr]) {
    GWASfunc(chr, index)
    print(index)
  }
}
```

**22. Filter for the same in IGAP set**

```
library(data.table)
GWAS_IGAPfilter <- function(chr, geneindex) {
  tryCatch({
    GWASname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_rawADGWAS/chr", chr, "/UKB_GWAS_chr", chr, "_gene", geneindex, ".rds")
    GWAS_result <- readRDS(GWASname)
    snpnamelist <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".snplist")
    snplist <- unlist(fread(snpnamelist, header = F))
    result <- GWAS_result[which(GWAS_result$rsID %in% snplist), ]
    GWASname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_IGAPset_ADGWAS/chr", chr, "/UKB_IGAPset_ADGWAS_chr", chr, "_gene", geneindex, ".rds")
    saveRDS(result, GWASname)
  }, error = function(e){print("error")})
}

gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)

for (chr in 1:22) {
  for (geneindex in 1:gene_n_v[chr]) {
    GWAS_IGAPfilter(chr, geneindex)
    print(geneindex)
  }
}
```

**23. Obtain 1000G SNPs by chr**
```
library(data.table)
KG_SNP_func <- function(chr) {
  chr_path <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/clump_0.5/chr", chr, "/", sep = "")
  all_files <- list.files(chr_path)
  matching_files <- grep(paste0("_chr", chr, "_0.5.clumps"), all_files, value = TRUE)
  
  clumped_ID <- c()
  for (i in 1:length(matching_files)) {
    setwd(chr_path)
    clumped <- fread(matching_files[i])
    clumped_rs <- unlist(clumped$ID)
    clumped_ID <- c(clumped_ID, clumped_rs)
    clumped_ID <- unique(clumped_ID)
  }
  outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/1000G_SNPlist_byChr/1000G_SNPlist_chr", chr, ".rds")
  saveRDS(clumped_ID, outname)
}
for (chr in 1:22) {
  KG_SNP_func(chr)
}
```

**24. Only intersect with 1000G**
```
library(data.table)
GWAS_1000Gfilter <- function(chr, geneindex) {
  tryCatch({
    GWASname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_rawADGWAS/chr", chr, "/UKB_GWAS_chr", chr, "_gene", geneindex, ".rds")
    GWAS_result <- readRDS(GWASname)
    snpnamelist <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/1000G_SNPlist_byChr/1000G_SNPlist_chr", chr, ".rds")
    snplist <- readRDS(snpnamelist)
    result <- GWAS_result[which(GWAS_result$rsID %in% snplist), ]
    GWASname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_1000Gset_ADGWAS/chr", chr, "/UKB_1000Gset_ADGWAS_chr", chr, "_gene", geneindex, ".rds")
    saveRDS(result, GWASname)
  }, error = function(e){print("error")})
}

gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)

for (chr in 1:22) {
  for (geneindex in 1:gene_n_v[chr]) {
    GWAS_1000Gfilter(chr, geneindex)
    print(geneindex)
  }
}
```

## Data check

**25. Check if snps are consistent in all sets**

Use *datacheck.R*

## Final data

**26. UKB IDP GWAS**

located in */gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5*, by chr, gene and IDP


**27. 1000G LD**

located in */gpfs/fs0/scratch/j/junpark/tianyu47/1000G/R_LD_0.5*, by chr and gene

**28. IGAP data**

located in */gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_0.5_bygene/*, by chr and gene

**29. IGAP snpset UKB GWAS**

located in */gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_IGAPset_ADGWAS*, by chr and gene

**29. UKB only version UKB GWAS**

located in */gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_1000Gset_ADGWAS*, by chr and gene
