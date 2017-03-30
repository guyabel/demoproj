library("tidyverse")

df0 <- read_csv("./data-raw/sweden1993_v0.csv")

df1 <- df0 %>%
  #work out survivorship rates and then project population
  mutate(sx_f = lead(Lx_f)/Lx_f,
         sx_m = lead(Lx_m)/Lx_m) %>%
  #surivorship for the penultimate age group
  mutate(sx_f = ifelse(test = row_number() == n() - 1,
                       yes = lead(Lx_f) / (lead(Lx_f) + Lx_f),
                       no = sx_f),
         sx_m = ifelse(test = row_number() == n() - 1,
                       yes = lead(Lx_m) / (lead(Lx_m) + Lx_m),
                       no = sx_m),
         #borrow the surivorship from the penultimate age group for projection
         #we could use lag() next, but sometimes we do know the last sx
         sx_f = ifelse(test = row_number() == n(),
                       yes = lag(sx_f),
                       no = sx_f),
         sx_m = ifelse(test = row_number() == n(),
                       yes = lag(sx_m),
                       no = sx_m)) %>%
  rename(fx = Fx)

sweden1993 <- df1
rm(df0, df1)
save(sweden1993, file = "data/sweden1993.rda")
