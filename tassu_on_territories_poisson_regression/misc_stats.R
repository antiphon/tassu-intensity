# stats

# inside 2011-2019
terri <- terr %>% 
  filter(start.tyear < 2020 & end.tyear >= 2011) %>%
  mutate(start.tyear = pmax(start.tyear, 2011),
         end.tyear  = pmin( end.tyear, 2020),
         duration.tyear = end.tyear - start.tyear)


terri %>% 
  #st_drop_geometry() %>% 
  group_by(type) %>% 
  summarise(total_time = sum(duration.tyear),
            n = n(),
            total_area = sum(st_area(geometry))/1000^2) 
