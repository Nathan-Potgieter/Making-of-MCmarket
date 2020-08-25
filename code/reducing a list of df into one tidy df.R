library(purrr)
1:100 %>% map(~rnorm(n = 50, mean = 0, sd = 1)) %>% reduce()


replicate(n = 100,
          expr = rnorm(n = 50, mean = 0, sd = 1),
          simplify = FALSE )


# how to use reduce
library(tidyverse)
collapsor <- function(X, ...) data.frame(value = rnorm(n = 50, mean = 0, sd = 1)) %>% dplyr::mutate(ID = X)
1:100 %>% map(~collapsor(X = ., n = 50, mean = 0, sd = 1)) %>% reduce(bind_rows)



