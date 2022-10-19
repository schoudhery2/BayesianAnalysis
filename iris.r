options(width=200)

data = read.table("iris.data",head=F,sep=',')
colnames(data) = c("sepal_len","sepal_wid","petal_len","petal_wid","species")

library(corrplot)
vals = data[,1:4]
corrplot(cor(vals)) # sepal_len and sepal_wid are least correlated


# just looking at predicting sepal_wid based on other 3 feats, independing of species

mod1 = lm(sepal_len ~ sepal_wid+petal_len+petal_wid,data=data)
pred1 = predict(mod1,data[,2:4])
plot(data$sepal_len,pred1)

#########################

for (sp in unique(data$species))
{
  cat(sprintf("mean(sepal_len) for %s = %s\n",sp,mean(data[data$species==sp,"sepal_len"])))
}

print(head(data))

temp = data
temp$versicolor = (data$species=="Iris-versicolor")
temp$setosa = (data$species=="Iris-setosa")
temp$virginica = (data$species=="Iris-virginica")


# dummy coding

coded = data
coded$versicolor = 0
coded[temp$versicolor==T,"versicolor"] = 1
coded$virginica = 0
coded[temp$virginica==T,"virginica"] = 1
coded$setosa = 0
coded[temp$setosa==T,"setosa"] = 1 # technically, I should choose one of these species as ref and set to 0


print(tail(coded))

# without species, there is an insignifcant relation between sepal_wid and sepal_len
mod2 = lm(sepal_len ~ sepal_wid,data=data)
print("mod2")
print(summary(mod2))

# R does dummy coding by default; one level is NA; R^2=0.72

mod3 = lm(sepal_len ~ sepal_wid+species,data=data)
print("mod3")
print(summary(mod3))

pred3 = predict(mod3,data)
plot(data$sepal_len,pred3)

# different intercept/offset for each species; R^2=0.99 (???)

mod3b = lm(sepal_len ~ 0+sepal_wid+species,data=data)
print("mod3b")
print(summary(mod3b))

pred3b = predict(mod3b,data)
plot(data$sepal_len,pred3b)

# make virginica the default factor level

data2 = data
data2$species = relevel(data2$species,ref="Iris-versicolor")
mod3c = lm(sepal_len ~ sepal_wid+species,data=data2)
print("mod3c")
print(summary(mod3c))


# using dummy coding - note how coeff for versicolor is NA

mod4 = lm(sepal_len ~ sepal_wid+setosa+virginica+versicolor,data=coded)
print("mod4")
print(summary(mod4))

# no independent intercept; coeffs are approximately the means above

mod5 = lm(sepal_len ~ 0+sepal_wid+setosa+virginica+versicolor,data=coded)
print("mod5")
print(summary(mod5))

# mod5 is same as mod3b

##########################

# contrast/effects/deviation coding, except I don't have a reference level (which should be excluded)

coded2 = data
coded2$versicolor = -1/3
coded2[temp$versicolor==T,"versicolor"] = 2/3
coded2$virginica = -1/3
coded2[temp$virginica==T,"virginica"] = 2/3
coded2$setosa = -1/3
coded2[temp$setosa==T,"setosa"] = 2/3

print(head(coded2))


mod6 = lm(sepal_len ~ 0+sepal_wid+setosa+virginica+versicolor,data=coded2)
print("mod6")
print(summary(mod6))

pred6 = predict(mod6,coded2)
plot(coded2$sepal_len,pred6)


# could try sum-coding, where reference level gets -1...
