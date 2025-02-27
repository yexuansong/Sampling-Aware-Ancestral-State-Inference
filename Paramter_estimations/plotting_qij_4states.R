r1 <- data1raw # no bias
r2 <- data2raw # state 2 2x oversampled
r3 <- data3raw # state 2 5x oversampled
r4 <- data4raw # state 2 10x oversampled
r5 <- data5raw # state 2 50x oversampled
r7 <- data7raw # state 2 100x oversampled, more time


m1 <- median(r1$x10)/0.1
m2 <- median(r2$x10)/0.1
m3 <- median(r3$x10)/0.1
m4 <- median(r4$x10)/0.1
m5 <- median(r5$x10)/0.1
m7 <- median(r7$x10)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  ggtitle("q41 estimate")+
  theme_minimal() +
  labs(x = "Sampling Ratio", y = "Median / true value")

m1 <- median(r1$x3)/0.1
m2 <- median(r2$x3)/0.1
m3 <- median(r3$x3)/0.1
m4 <- median(r4$x3)/0.1
m5 <- median(r5$x3)/0.1
m7 <- median(r7$x3)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q14 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x12)/0.1
m2 <- median(r2$x12)/0.1
m3 <- median(r3$x12)/0.1
m4 <- median(r4$x12)/0.1
m5 <- median(r5$x12)/0.1
m7 <- median(r7$x12)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))


ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q14 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")

m1 <- median(r1$x1)/0.1
m2 <- median(r2$x1)/0.1
m3 <- median(r3$x1)/0.1
m4 <- median(r4$x1)/0.1
m5 <- median(r5$x1)/0.1
m7 <- median(r7$x1)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q12_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q12 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")

m1 <- median(r1$x2)/0.1
m2 <- median(r2$x2)/0.1
m3 <- median(r3$x2)/0.1
m4 <- median(r4$x2)/0.1
m5 <- median(r5$x2)/0.1
m7 <- median(r7$x2)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q13_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q13 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")

m1 <- median(r1$x3)/0.1
m2 <- median(r2$x3)/0.1
m3 <- median(r3$x3)/0.1
m4 <- median(r4$x3)/0.1
m5 <- median(r5$x3)/0.1
m7 <- median(r7$x3)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q14_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q14 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x4)/0.1
m2 <- median(r2$x4)/0.1
m3 <- median(r3$x4)/0.1
m4 <- median(r4$x4)/0.1
m5 <- median(r5$x4)/0.1
m7 <- median(r7$x4)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q21_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q21 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")

m1 <- median(r1$x5)/0.1
m2 <- median(r2$x5)/0.1
m3 <- median(r3$x5)/0.1
m4 <- median(r4$x5)/0.1
m5 <- median(r5$x5)/0.1
m7 <- median(r7$x5)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q23_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q23 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x6)/0.1
m2 <- median(r2$x6)/0.1
m3 <- median(r3$x6)/0.1
m4 <- median(r4$x6)/0.1
m5 <- median(r5$x6)/0.1
m7 <- median(r7$x6)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q24_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q24 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x7)/0.1
m2 <- median(r2$x7)/0.1
m3 <- median(r3$x7)/0.1
m4 <- median(r4$x7)/0.1
m5 <- median(r5$x7)/0.1
m7 <- median(r7$x7)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q31_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q31 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x8)/0.1
m2 <- median(r2$x8)/0.1
m3 <- median(r3$x8)/0.1
m4 <- median(r4$x8)/0.1
m5 <- median(r5$x8)/0.1
m7 <- median(r7$x8)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q32_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q32 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x9)/0.1
m2 <- median(r2$x9)/0.1
m3 <- median(r3$x9)/0.1
m4 <- median(r4$x9)/0.1
m5 <- median(r5$x9)/0.1
m7 <- median(r7$x9)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q34_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q34 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x10)/0.1
m2 <- median(r2$x10)/0.1
m3 <- median(r3$x10)/0.1
m4 <- median(r4$x10)/0.1
m5 <- median(r5$x10)/0.1
m7 <- median(r7$x10)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q41_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q41 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x11)/0.1
m2 <- median(r2$x11)/0.1
m3 <- median(r3$x11)/0.1
m4 <- median(r4$x11)/0.1
m5 <- median(r5$x11)/0.1
m7 <- median(r7$x11)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q42_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q42 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


m1 <- median(r1$x12)/0.1
m2 <- median(r2$x12)/0.1
m3 <- median(r3$x12)/0.1
m4 <- median(r4$x12)/0.1
m5 <- median(r5$x12)/0.1
m7 <- median(r7$x12)/0.1

plot_data <- data.frame(ratio = c(1,2,5,10,50,100), median = c(m1,m2,m3,m4,m5,m7))

q43_plot<-ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 50, 100)) +
  theme_minimal() + ggtitle("q43 estimate")+
  labs(x = "Sampling Ratio", y = "Median / true value")


pp <- ggarrange(q12_plot,q13_plot,q14_plot,q21_plot,q23_plot,q24_plot,q31_plot,
                q32_plot,q34_plot,q41_plot,q42_plot,q43_plot, nrow=4,ncol=3)

