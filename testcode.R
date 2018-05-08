rm(list = ls())

install_local(wd())


require(SwitchPack)

sim.df <- simStudy()

x <- rpsft.input(Surv(os.t, os.e) ~ I(x.trt==1),
                 Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
                 AdminCensTime = t.censor,
                 data = sim.df)

plot(x)
y <- RPSFT(x)

y <- RPSFT.latent(0.4, x)


require(dplyr)
require(Hmisc)
require(GGally)


survfit(Surv(pfs.t, event = pfs.e) ~ x.trt, data = sim.df) %>%
  plot()



survfit(Surv(os.t, event = os.e) ~ x.trt, data = sim.df) %>%
  plot()








latent <- RPSFT.latent(psi = 0.4, rpsft.input = x)




survfit(Surv(os.t, event = os.e) ~ x.trt, data = sim.df) %>%
  ggsurv(back.white = TRUE)

survfit(Surv(os.t, event = os.e) ~ x.trt, data = sim.df) %>%
  ggsurv(back.white = TRUE) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw()


RPSFT.cox()

require(reshape2)

v.df <- rbind(data_frame(value = c(100,50,25,10), cost = c(113,58,31,15), shop = "Media Markt"),
              data_frame(value = c(100,50,20,10), cost = c(113,59,26,15), shop = "Manor"))

evalCV <- function(x, df = v.df){
  df %>%
  mutate(n.sel = c(x,x), tvalue = value*n.sel, tcost = cost*n.sel) %>%
  group_by(shop) %>%
  summarise(cost = sum(tcost), value = sum(tvalue))
}


combs <- expand.grid(0:2,0:3,0:4,0:8) %>%
  mutate(cid = 1:540)

rc.df <- evalCV(as.numeric(combs[1,1:4])) %>%
  mutate(cid = 1)

for (i in 2:540){
  this.df <- evalCV(as.numeric(combs[i,1:4])) %>%
    mutate(cid = i)
  rc.df <- rbind(rc.df, this.df)
}

check.df <- rc.df %>%
  filter(cost < 194) %>%
  group_by(shop) %>%
  arrange(desc(value),cost)

require(ggplot2)

ggplot(rc.df, aes(x=cost, y = value, color = shop, shape = shop)) +
  geom_point()


cid.sel <- rc.df %>%
  filter(value == 160, cost <188) %>%
  transmute(cid) %>%
  unlist() %>%
  as.numeric()

combs %>%
  filter(cid %in% cid.sel)

test.df <- data_frame(trt = rep(c("Pla","Act1","Act2"), each = 3),Attr=rep(c("Attr1","Attr2","Attr3"),3), value = runif(9))
n.attr <- length(unique(test.df$Attr))
angle.off <- 2*pi / (n.attr+1)
df <- mutate(test.df, x = angle.off*c(0:2))

ggplot(df, aes(x=x, y=value, fill = trt))+
  geom_line(color = "black") +
  geom_polygon(alpha =  0.5) +
  geom_point() +
  coord_polar(theta = "x")

evalCV(x1)
