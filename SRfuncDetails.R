# Smsy
SR <- function(ss,pp,cc) ss/(1/pp^ww + (ss/cc)^ww)^(1/ww)
Smsy <- function(pp,cc,ww=3) cc/pp*(pp^ww/(ww+1)-1)
Umsy <- function(pp,ww=3) 1-pp^(-ww/(ww+1))
pp <- seq(1,10,length.out=100)

plot(pp,Umsy(pp,ww=1),type="l",lwd=2,ylim=c(0,1))
for(ww in 2:100){
  lines(pp,Umsy(pp,ww))
}
lines(pp,Umsy(pp,ww=1000),col="blue",lwd=2)

# d/ds(s*((1/(p^w+(s/c)^w)^(1/w))-1)) in wolfram alpha (notice p should be 1/p)
# an alternate solution is:
# p^w *((s/c)^w + p^w)^(-(w + 1)/w) - 1
# setting to zero and solving gives
# Solve[-1 + p^w/(p^w + (s/c)^w)^((1 + w)/w) == 0, {s}]
# s = c ((p^(-w))^(-w/(w + 1)) - p^w)^(1/w)
# substituting 1/p for p gives
# s = c ((p^(-w*w/(w + 1)) - p^(-w))^(1/w) # S at MSY
# putting this into (f(s)-s)/f(s) = 1-s/f(s) produces
# 1-p^(-w/(w+1)) # U

# steepness is 
# p * (1-0.2^w+(0.2*p)^w)^(-1/w)
# 
