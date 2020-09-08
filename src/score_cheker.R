library("tidyverse")
setwd("~/work/irabu")

data <- read_delim(delim = " ","./result/negative_control.csv",col_names = FALSE)
offline <- read_csv("../negative_control/result/scores_mock.csv",col_names= FALSE)
data <- data %>% rename(score = X2) %>% mutate(type = "irabu") %>% filter(score != 0) %>% select(type,score) %>% bind_rows(offline %>% rename(score = X1) %>% mutate(type = "offline"))

summary <- data %>% nest(-type) %>% mutate(data = map(data,function(df) df %>% summarize(mean = mean(score),sd = sd(score),max = max(score),min = min(score)))) %>% unnest()


