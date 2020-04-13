
B =1000000
Time = 100
X = matrix(100, B, Time)
X2 = matrix(100, B, Time)
Y = matrix(100, B, Time)
Y2 = matrix(100, B, Time)
R0  = 1
m = 1.5
for (b in 1:B){
  for (t in 1:(Time-1)){
    X[b, t+1] = max(c(0,(rnorm(1, 2 *sqrt(3/8 + R0 * X[b, t]),1)/2)^2 - 3/8))
    X2[b, t+1] = max(c(0,(rnorm(1, 2 *sqrt(3/8 + rgamma(1, R0 * 10, 10) * X2[b, t]),1)/2)^2 - 3/8))
    Y[b, t+1] = max(c(0,(rnorm(1, 2 *sqrt(3/8 + m * R0 * Y[b, t]),1)/2)^2 - 3/8))
    Y2[b, t+1] = max(c(0,(rnorm(1, 2 *sqrt(3/8 + rgamma(1, m * R0 * 10, 10) * Y2[b, t]),1)/2)^2 - 3/8))
    
      }
}

X  = X[1:425500,]
X2  = X2[1:425500,]
Y  = Y[1:425500,]
Y2  = Y2[1:425500,]

stopping_time <- apply(X, 1, FUN= function(x)((cumsum(x)[Time]> 5000 * 100)))
print(c(100 * mean(apply(X, 1, FUN= function(x)((cumsum(x)[Time]> 5000 * 100)))),
        100 * mean(apply(X2, 1, FUN= function(x)((cumsum(x)[Time]> 5000 * 100)))),
        100 * mean(apply(Y, 1, FUN= function(x)((cumsum(x)[Time]> 5000 * 100)))),
        100 *  mean(apply(Y2, 1, FUN= function(x)((cumsum(x)[Time]> 5000 * 100))))))

print(c(mean(apply(t(apply(X2,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1])), na.rm=TRUE) ,
        mean(apply(t(apply(Y,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1]))),
        mean(apply(t(apply(Y2,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1])))
        )
)
stopping_time <- apply(t(apply(X2,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1]))
stopping_time2 <- apply(t(apply(Y,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1]))
stopping_time3 <- apply(t(apply(Y2,1,cumsum)), 1, FUN= function(x)(which(x > 5000 * 100)[1]))
print(c(mean(stopping_time, na.rm=TRUE), mean(stopping_time2),mean(stopping_time3)))      

stopping_time = data.frame(y= stopping_time[!is.na(stopping_time)],
                           type = rep("Random R=Gamma(1,1)",length(stopping_time[!is.na(stopping_time)]))
                           )
stopping_time = rbind(stopping_time,
                      data.frame(y = stopping_time2[!is.na(stopping_time2)],
                                 type = rep("Constant R=1.5",length(stopping_time2[!is.na(stopping_time2)]))
                      ),
                      data.frame(y = stopping_time3[!is.na(stopping_time3)],
                                 type = rep("Random R=Gamma(1.5,1)",length(stopping_time3[!is.na(stopping_time3)]))
                      )
                      )
ggplot(stopping_time)+
  geom_bar(aes(x=y, y=..y../sum(..y..), fill=type)) +theme_bw()
df = data.frame(matrix(0, 100,11))
df2 = data.frame(matrix(0, 100,11))
df3 = data.frame(matrix(0, 100,11))
df4 = data.frame(matrix(0, 100,11))
quants = c(0.001,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.999)
for (q in 1:11){
  df[,q] = apply(t(apply(X,1,cumsum)), 2, FUN=function(x){quantile(x,quants[q])})
  df2[,q] = apply(t(apply(X2,1,cumsum)), 2, FUN=function(x){quantile(x , quants[q])})
  df3[,q] = apply(t(apply(Y,1,cumsum)), 2, FUN=function(x){quantile(x, quants[q])})
  df4[,q] = apply(t(apply(Y2,1,cumsum)), 2, FUN=function(x){quantile(x, quants[q])})
}
df[,12] = apply(t(apply(X,1,cumsum)), 2, mean)
df2[,12] = apply(t(apply(X2,1,cumsum)), 2, mean)
df3[,12] = apply(t(apply(Y,1,cumsum)), 2, mean)
df4[,12] = apply(t(apply(Y2,1,cumsum)), 2,mean)
df$x =1:100
df2$x =1:100
df3$x =1:100
df4$x =1:100
df$type = "Constant R=1"
df2$type = "Random R=Gamma(1,1)"
df3$type = "Constant R=1.5"
df4$type = "Random R=Gamma(1.5,1)"

df = rbind(df, df2,df3,df4)
for (q in 1:11){
  colnames(df)[q] = paste("quantile",toString(100 * quants[q]))
}
colnames(df)[12] ="mean"
df$type = as.factor(df$type)
ggplot(df %>%filter(type %in% c("Constant R=1","Random R=Gamma(1,1)")))+
  geom_point(aes(x=x, y = `mean` ),col="black", shape=18, size=3)+
  geom_errorbar(aes(x=x,ymin = `quantile 1`, ymax = `quantile 99`, col=type),position=pd)+
  geom_point(aes(x=x, y = `quantile 50`, col=type ), shape=20, size=2)+
  theme_bw()+scale_y_continuous(trans='log10')+
  theme(legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=16),
                                                                axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),
                                                                axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
ggsave("~/Desktop/mean_simr1.pdf")
points(log(apply(t(apply(X2,1,cumsum)), 2, mean)), col='red', pch=3, m=1)

plot( a)