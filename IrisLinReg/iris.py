import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf

data_df = pd.read_csv("./iris.data", header=None)
data_df.columns = ["sepal_len","sepal_wid","petal_len","petal_wid","species"]


vals = data_df.drop(columns=["species"])
corr = vals.corr()

# Generate a mask for the upper triangle
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True
# sns.heatmap(corr, mask=mask, cmap=plt.cm.PuOr,square=True)
# plt.show()


# just looking at predicting sepal_wid based on other 3 feats, independing of species

mod1 = smf.ols(formula='sepal_len ~ sepal_wid + petal_len + petal_wid', data=data_df).fit()
pred1 = mod1.predict(data_df[['sepal_wid','petal_len','petal_wid']])
plt.scatter(data_df["sepal_len"],pred1, s=4, label = "mod1")

#plt.show()

#########################

for sp in set(data_df["species"]):
  print("mean(sepal_len) for %s = %s"%(sp,np.mean(data_df[data_df["species"]==sp]["sepal_len"])))

temp = data_df
temp["versicolor"] = data_df["species"]=="Iris-versicolor"
temp["setosa"] = data_df["species"]=="Iris-setosa"
temp["virginica"] = data_df["species"]=="Iris-virginica"

# dummy coding

coded = data_df
coded["versicolor"] = 0
coded[temp["versicolor"]==True]["versicolor"] = 1
coded["virginica"] = 0
coded[temp["virginica"]==True]["virginica"] = 1
coded["setosa"] = 0
coded[temp["setosa"]==True]["setosa"] = 1 # technically, I should choose one of these species as ref and set to 0

print(coded)
# without species, there is an insignifcant relation between sepal_wid and sepal_len
mod2 = smf.ols(formula='sepal_len ~ sepal_wid',data=data_df).fit()
print("mod2")
print(mod2.summary())

# R does dummy coding by default; one level is NA; R^2=0.72

mod3 = smf.ols(formula='sepal_len ~ sepal_wid+species',data=data_df).fit()
print("mod3")
print(mod3.summary())

pred3 = mod3.predict(data_df)
plt.scatter(data_df["sepal_len"],pred3, s=4, label = "mod3")

# different intercept/offset for each species; R^2=0.99 (???)

mod3b = smf.ols(formula='sepal_len ~ 0+sepal_wid+species',data=data_df).fit()
print("mod3b")
print(mod3b.summary())

pred3b = mod3b.predict(data_df)
plt.scatter(data_df["sepal_len"],pred3b, s=4, label = "mod3b")



# make virginica the default factor level

data2 = data_df
mod3c = smf.ols(formula='sepal_len ~ sepal_wid + C(species,Treatment(reference="Iris-versicolor"))',data=data2).fit()
print("mod3c")
print(mod3c.summary())

# using dummy coding - note how coeff for versicolor is NA

mod4 = smf.ols(formula='sepal_len ~ sepal_wid+setosa+virginica+versicolor',data=coded).fit()
print("mod4")
print(mod4.summary())

# no independent intercept; coeffs are approximately the means above

mod5 = smf.ols(formula='sepal_len ~ 0+sepal_wid+setosa+virginica+versicolor',data=coded).fit()
print("mod5")
print(mod5.summary())

# mod5 is same as mod3b

##########################

# contrast/effects/deviation coding, except I don't have a reference level (which should be excluded)

coded2 = data_df
coded2["versicolor"] = -1/3
coded2[temp["versicolor"]==True]["versicolor"] = 2/3
coded2["virginica"] = -1/3
coded2[temp["virginica"]==True]["virginica"] = 2/3
coded2["setosa"] = -1/3
coded2[temp["setosa"]==True]["setosa"] = 2/3

print(coded2)


mod6 = smf.ols(formula='sepal_len ~ 0+sepal_wid+setosa+virginica+versicolor',data=coded2).fit()
print("mod6")
print(mod6.summary())

pred6 = mod6.predict(coded2)
plt.scatter(coded2["sepal_len"],pred6, label="mod6", s=4)

# could try sum-coding, where reference level gets -1...

plt.legend()
plt.title("Linear Regression Models Testing Various Coding")
plt.xlabel("Actual Sepal Length")
plt.ylabel("Predicted Sepal Length")
plt.show()