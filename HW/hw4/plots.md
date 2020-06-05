---
title: "Extra Plots"
author: "Ethan Shen"
date: "6/4/2020"
output: 
  html_document:
    keep_md: TRUE
---




```r
nc112 = st_read('data/nc_districts112.shp', quiet=TRUE, stringsAsFactors=TRUE)
nc114 = st_read('data/nc_districts114.shp', quiet=TRUE, stringsAsFactors=TRUE)
nc=st_read("nc_counties.gpkg", quiet=TRUE)
```


# 3C.1


```r
ggplot() + 
  geom_sf(data=nc112, aes(fill=DISTRICT)) + 
  geom_sf(data=nc, fill=NA) + 
  labs(title = "NC112") + 
  guides(fill = guide_legend(ncol=2))
```

![](plots_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
ggplot() + 
  geom_sf(data=nc114, aes(fill=DISTRICT)) + 
  geom_sf(data=nc, fill=NA) + 
  labs(title = "NC114") + 
  guides(fill = guide_legend(ncol=2))
```

![](plots_files/figure-html/unnamed-chunk-2-2.png)<!-- -->



```r
nc112_intersection = st_intersection(nc, nc112  %>%
                                       st_transform(st_crs(nc))%>% filter(DISTRICT==11))

nc114_intersection = st_intersection(nc, nc114  %>%
                                       st_transform(st_crs(nc))%>% filter(DISTRICT==11))
ggplot()  + 
  geom_sf(data=nc112_intersection,
          aes(fill=COUNTY)) +
  geom_sf(data = nc112 %>% filter(DISTRICT==11), fill=NA)
```

![](plots_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ggplot() + 
  geom_sf(data=nc114_intersection,
          aes(fill=COUNTY))+
  geom_sf(data = nc114 %>% filter(DISTRICT==11) , fill=NA)
```

![](plots_files/figure-html/unnamed-chunk-3-2.png)<!-- -->
