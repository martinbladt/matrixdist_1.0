
x <- new("ph", name = "my_ph")
x
x@name
x@pars
x@fit

coef(x)

y <- sim(x)

f <- fit(x,y)
f@fit

