library(faraway)
data(seatpos)

# "the pita breads are not the shadows of the footballs"

wt_cat <- cut(seatpos$Weight, include.lowest=TRUE, breaks=summary(seatpos$Weight)[-4])

plot_ly (type = "scatter3d",
         x = seatpos$Arm, 
         y = seatpos$Thigh,
         z = seatpos$hipcenter,
         mode = "markers", 
         color = wt_cat) %>% 
  layout(scene = list(xaxis = list(title = 'Arm'),
                      yaxis = list(title = 'Thigh'),
                      zaxis = list(title = 'hip')),
         title = 'by weight categories')


library(readr)
posterior <- read_csv("posterior_summary/thinned_MCMC.csv") # x1 is the row name


select(posterior, starts_with("ind"),X1) %>%
  gather(key, value, -X1) %>%
  group_by(key) %>%
  summarise(non_zero_freq =  mean(value))

# ind from 1 to 20, alpha (?), betaT 1 to 20

with(posterior, plot_ly (type = "scatter3d",
                         x = `betaT[9]`*`Ind[9]`, 
                         y = `betaT[10]`*`Ind[10]`,
                         z = `betaT[20]`*`Ind[20]`,
                         mode="markers",
                         opacity = 0.01)) %>%
  layout(scene = list(xaxis = list(title = '9'),
                      yaxis = list(title = '10'),
                      zaxis = list(title = '20'),
                      title = 'points only'))


## need an object that stores all things in one plotly object
## select three variables
### for each draw, how many are non-zero and in which orders

is.zero <- function(x){
  
  return(sum((x != 0) * c(4,2,1)))
  
}



to_plot <- transmute(posterior, 
                     x = `betaT[5]`*`Ind[5]`, 
                     y = `betaT[9]`*`Ind[9]`,
                     z = `betaT[15]`*`Ind[15]`) %>%
  mutate(sample =1:nrow(.)) %>%
  split(.$sample) %>%
  map_df(~mutate(.x, indicator = is.zero(.x[,1:3])))

to_plot %<>% mutate(bread = c("Donut hole",
                              "Grissini 1",
                              "Grissini 2",
                              "Pita 1",
                              "Grissini 3",
                              "Pita 2",
                              "Pita 3",
                              "Dinner roll")[indicator + 1])

with(to_plot, plot_ly (type = "scatter3d",
                         x = x, 
                         y = y,
                         z = z,
                         mode="markers",
                       color=bread,
                       opacity = 0.25)) %>%
  layout(scene = list(xaxis = list(title = '5'),
                      yaxis = list(title = '9'),
                      zaxis = list(title = '15'),
                      title = 'points only'))


windows()
select(posterior, starts_with('beta')) %>%
  mutate(idx =1:nrow(.)) %>%
  gather(key, value, -idx) %>%
  ggplot(data=., aes(x=idx, y=value)) +
  geom_line() +
  facet_wrap(~ key)


select(posterior, starts_with('Ind')) %>%
  mutate(idx =1:nrow(.)) %>%
  gather(key, value, -idx) %>%
  group_by(key) %>%
  summarise(Inclusion = mean(value)) %>%
  mutate(realkey = parse_number(key)) %>%
  ggplot(data=., aes(x=realkey, y=Inclusion)) +
  geom_col()


