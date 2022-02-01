# Brian Tinsley
# Week 2 Homework

#exercise 1.1
#1. Identify the rows in which there was no reported station in the station column. (Hint: use is.na()).
which(is.na(attenu$station))
#or
attenu[is.na(attenu$station),]
#2. Create a copy of attenu called attenu_cleaned, where rows that are missing station information are
#not included.
attenu_cleaned <- attenu[!is.na(attenu$station),] 
attenu_cleaned
#3. Print the first 6 rows of attenu_cleaned (using head()) and its dimensions (using dim()).
head(attenu_cleaned)
dim(attenu_cleaned)


#exercise 1.2
#1. Make a copy of the Theoph data set called Theoph_2 (again, the data set should be included in R).
Theoph_2 <- Theoph
#2. Identify what the median treatment dose is. (Use str() to look see which column has dose information).
med <- median(Theoph_2$Dose)
#3. Then, using ifelse(), add a new column called Dose_Class to Theoph_2, where the dose is classified
#as high if it is above or equal to the median dose, and 'low otherwise. (Use the $ to add the column).
Theoph_2$Dose_Class <- ifelse(Theoph_2$Dose >= med, "high", "low")
#4. As before, print the first 6 rows and its dimensions.
head(Theoph_2)
dim(Theoph_2)


#exercise 1.3
#1. Read in the Starbucks nutrition data using read.csv() and store it in a variable named 
#   starbucks. The file is called starbucks.csv. As with using the command line, make sure that 
#   you are in the correct directory by using getwd() and setwd()!
starbucks <- read.csv('starbucks.csv')
#2. Now, let’s clean up our data by filtering out the drinks without any data. Follow these steps:
#  a. Determine which values contain NA using is.na().
starb_na <- is.na(starbucks)
#  b. Notice that we get a data frame of TRUE and FALSE data. How can we identify which 
#     rows are full of NA values? Create a boolean vector, is_row_empty to identify these rows. 
#     (Hint: remember that TRUE is 1 and FALSE is 0. What function sums across rows? 
#     How many NAs in a row do we need to say that there is no nutritional information for that drink?)
#
is_row_empty <- rowSums(starb_na) #a 6 means the row is empty
#  c. Use nrow() or dim() to identify the number of rows in starbucks, and verify that is_row_empty 
#     is the same length as the number of rows. (Make sure length(is_row_empty) == nrow(starbucks)).
nrow(starbucks) == length(is_row_empty)

#  d. With the boolean vector, create starbucks_cleaned, where the drinks without nutritional data 
#     are removed. (Hint: does the boolean vector select rows or columns? How do we select rows from a data frame?)
starbucks_cleaned <- starbucks[!(is_row_empty == 6),]
#3. Now that we’ve cleaned everything up, let’s examine the data some more. Create a plot of Calories (y- axis) vs. Carb. 
#   What do we notice about the plot? Make sure the rename the axes! (Note: carbohydrates are in grams; be sure to specify 
#   that in the name.)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab="Carbs (grams)",ylab="Calories")
#the more carbs there are, the more calories there are
#4. Let’s look at the drink with the highest amount of calories. Identify what is the highest number of calories 
#   by using the max() function. Then, use boolean indexing to identify which row has this drink with the 
#   highest number of calories (you’ll get a 1 x 7 data frame). Find out the name of the drink by using either the 
#   square bracket [] or $ syntax.
starbucks_cleaned[starbucks_cleaned$Calories == max(starbucks_cleaned$Calories),]
#5. Let’s highlight this point in the plot from part (4). Add a new column, is_highest_fat, which contains TRUE 
#   if the drink has the highest amount of fat, and FALSE otherwise. Then, reusing the code from (3), color the 
#   plot points by this column.
starbucks_cleaned$is_highest_fat <- ifelse(starbucks_cleaned$Calories == max(starbucks_cleaned$Calories), TRUE, FALSE)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab="Carbs (grams)",ylab="Calories",col=factor(starbucks_cleaned$is_highest_fat))
#6. BONUS: now, color all points by a gradient based on their fat content. (Google how to do so! It’s a lot easier in ggplot.)


#exercise 1.4
#1. Read the data in using the read.csv() function.
batting <- read.csv('Batting.csv')
#2. Identify the number of players who scored three or more home runs in a given year. The number of
#   homeruns in in the HR column. (note: players are duplicated, but for the purposes of this question just
#   pretend that everyone is unique).
sum(batting$HR >= 3)
#3. Plot the number of homeruns vs. year. (Hint: write a function using the plot() function to not repeat
#   doing this analysis, as you’ll see in the next steps).
plot(batting$yearID, batting$HR)
#4. Create a new data frame, containing players from the LA Angels. This information is in the teamID
#   column, and the LA Angels code is “LAA”. Then, repeat this analysis, but this time restricting your
#   data to the LA Angels players.
la_angels <- batting[batting$teamID == "LAA",]
plot(la_angels$yearID, la_angels$HR)
#5. Repeat step 4, except subset the original batting data to include players from either “ATL” or “PIT”
#   (so, this step should be one data frame). Then, in the scatter plot, color the number of homeruns by the team.
atl_or_pit <- batting[batting$teamID == "PIT" | batting$teamID == "ATL",]
plot(atl_or_pit$yearID, atl_or_pit$HR, col=factor(atl_or_pit$teamID))


#exercise 1.5
easy_plot <- function(x,y,color_data){
#1. Identify the median data value of color_data. Then, using ifelse() create a new vector, levels, where values
#   less than the median are "low" and values higher than the median are "high". Then convert color_levels into factor 
#   by casting it as a factor as follows (you need to do this for coloring): levels = factor(levels). Print this median.
  medi <- median(color_data)
  levels <- ifelse(color_data < medi, "low", "high")
  levels = factor(levels)
  print(medi)
#2. Using plot() or ggplot(), plot x and y and color by color_levels. Change the points to filled-in circles using 
#   the option pch=20.
  plot(x,y,col=levels,pch=20)
#3. Finally, print the correlation between x and y (basically how related they are), by using the function cor.test(x, y). 
#   You should see the correlation and a p-value. Note: you won’t be returning anything in this function.
  print(cor.test(x,y))
}
#4. Test your function by passing the same x twice: to x and to color_data. How are the points colored?
easy_plot(starbucks_cleaned$Calories, starbucks_cleaned$Calories, starbucks_cleaned$Calories)
#points are colored by being less or greater than the median

#Now, use easy_plot() to examine the cleaned Starbucks and batting data sets (pick any three variables to
#plot – make sure they’re numeric). Do you see significant correlations between the variables you picked?
easy_plot(starbucks_cleaned$Protein, starbucks_cleaned$Fat, starbucks_cleaned$Calories)




#exercise 2.1
#What does the data set describe? How many observations does it contain (observations means data points)? 
#What features does the data set contain (features means the number of variables per the data point)?
head(iris)
?iris
#this data set gives measurements of sepal length and width in flowers from different species

#exercise 2.2 
#Which features are continuous variables? (Continuous variables are always numeric, but the reverse is not 
#necessarily true, especially when an integer is used encode a category. For example, it’s common for 1 to 
#represent category 1, 2 for category 2, etc.) Which are categorical? (Factors are categorical variables, 
#but categorical variables are not necessarily labeled as factors such in the data set). Can you tell me 
#the R data type for each column?
colnames(iris)
#species is categorical, sepal length/width, and petal length/width are all continuous doubles

#exercise 2.3
#Examine how each of the continuous variables look, making a histogram for each of them. 
#Do you notice anything interesting?
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)
#petal lengths/widths are more split while sepal length and width are more normally distrubuted

#exercise 2.4
#Let’s say we want to split the data set into “narrow-sepaled” plants and “wide-sepaled” plants. 
#To do so, let’s want to add a new column to the iris data set containing this data and perform 
#some more analysis. Follow these steps:
#   1. Find the mean sepal width in the data frame, and store it in a variable. Also, make a copy of iris called iris_copy.
mean_sepal <- mean(iris$Sepal.Width)
iris_copy <- iris
#   2. Create a new vector using the ifelse() function to fill out values, comparing Sepal.Width to the mean value.
sepal_vec <- ifelse(iris$Sepal.Width < mean_sepal, "low","high")
#   3. Create a new column in the iris data set (using the dollar sign) and assign the new vector to that column.
iris_copy$width_bool <- sepal_vec
#   4. Create a boxplot plotting the Sepal.Width based on the new column (narrow vs. wide). To plot a boxplot by a variable, 
#you need to use the syntax y ~ x instead of x = x, y = y as with plot(), where x and y are variables.
boxplot(iris_copy$Sepal.Width ~ iris_copy$width_bool)

#exercise 2.5
#Examine the following plots. The x and y axes of these plots are specified on the bottom and left, respectively
#Based on these plots, which species looks the most unique out of the three? Which species look the most similar?
#the species Setosa is the most unique while the other two are relatively similar
pairs(iris[1:4],pch = 21, bg = c("red", "green", "blue")[unclass(iris$Species)])
