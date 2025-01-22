#usdpt of education data which is supposed to be institution specific
#"PS" is privacy surpessed 

# look at schools like dartmouth vs schools like yale/harvard
# is there any measureable differnce bewtween two very similar schools
# picking 2-3 schools that are close to eachother on the us national rankings -- look at duke or stanford or somthin


library(tidyverse)
library(readr)
library(ggplot2)

# schools used ------------------------------------------------------------

picked_dartmouth<- c(
  "Johns Hopkins University", 
  "Rice University", 
  "Bucknell University", 
  "Wesleyan University", 
  "Bowdoin College", 
  "Cooper Union for the Advancement of Science and Art", 
  "Colgate University", 
  "Tufts University", 
  "College of William and Mary", 
  "Harvey Mudd College", 
  "University of Notre Dame", 
  "Wake Forest University", 
  "Alabama A and M University", 
  "Georgetown University", 
  "Barnard College", 
  "Yale University", 
  "University of Pennsylvania", 
  "Vanderbilt University", 
  "Stanford University", 
  "Harvard University", 
  "Brown University"
)


dart_picks <- c("Yale University",
                "Princeton University", 
                "Harvard University", 
                "Cornell University",
                "Stanford University", 
                "Massachusetts Institute of Technology",
                "Vanderbilt University",
                "University of Pennsylvania", 
                "University of Chicago", 
                "Columbia College",
                "Dartmouth College")


dartmouth_association<- c(picked_dartmouth, dart_picks)
dartmouth_association <- unique(dartmouth_association)
print(dartmouth_association)


# data sources ------------------------------------------------------------
dart_association <- c(picked_dartmouth, dart_picks)

most_recent <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/Most-Recent-Cohorts-Field-of-Study.csv")
m <- most_recent %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>%  # Convert the column to numeric
  filter(!is.na(EARN_MDN_HI_2YR))   # Filter out rows where the column is NA
  # filter(dart_association)


fod1<- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1617_1718_PP.csv")
fd1 <- fod1 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>%  
  filter(!is.na(EARN_MDN_HI_2YR))
  # filter(dart_association)
  # filter(picked_dartmouth)

fod2 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1819_1920_PP.csv")
fd2 <- fod2 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
  filter(!is.na(EARN_MDN_HI_2YR)) 

fod3 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1718_1819_PP.csv")
fd3 <- fod3 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
  filter(!is.na(EARN_MDN_HI_2YR))  

fod4 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1920_2021_PP.csv")
fd4 <- fod4 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
  filter(!is.na(EARN_MDN_HI_2YR))  

fod5 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1516_1617_PP.csv")
fd5 <- fod5 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
  filter(!is.na(EARN_MDN_HI_2YR))  

fod6 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/FieldOfStudyData1415_1516_PP.csv")
fd6 <- fod6 %>%
  mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
  filter(!is.na(EARN_MDN_HI_2YR))  

#mo97 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED1996_97_PP.csv")

# m97 <- mo97 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo98 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED1997_98_PP.csv")
# m98 <- mo98 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo99 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED1998_99_PP.csv")
# m99 <- mo99 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo00 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED1999_00_PP.csv")
# m00 <- mo00 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo01 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2000_01_PP.csv")
# m01 <- mo01 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo02 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2001_02_PP.csv")
# m02 <- mo02 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo03 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2002_03_PP.csv")
# m03 <- mo03 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo04 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2003_04_PP.csv")
# m04 <- mo04 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo05 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2004_05_PP.csv")
# m05 <- mo05 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo06 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2005_06_PP.csv")
# m06 <- mo06 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo07 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2006_07_PP.csv")
# m07 <- mo07 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo08 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2008_09_PP.csv")
# m08 <- mo08 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo09 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2008_09_PP.csv")
# m09 <- mo09 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo10 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2009_10_PP.csv")
# m10 <- mo10 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR))  
# 
# mo11 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2010_11_PP.csv")
# m11 <- mo11 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo12 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2011_12_PP.csv")
# m12 <- mo12 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo13 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2012_13_PP.csv")
# m13 <- mo13 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo14 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2013_14_PP.csv")
# m14 <- mo14 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo15 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2014_15_PP.csv")
# m15 <- mo15 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo16 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2015_16_PP.csv")
# m16 <- mo16 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo17 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2016_17_PP.csv")
# m17 <- mo17 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo18 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2017_18_PP.csv")
# m18 <- mo18 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo19 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2018_19_PP.csv")
# m19 <- mo19 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo20 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2019_20_PP.csv")
# m20 <- mo20 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo21 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2020_21_PP.csv")
# m21 <- mo21 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 
# 
# mo22 <- read_csv("us_dpt_ed_college_scorecard/College_Scorecard_Raw_Data_09242024/MERGED2021_22_PP.csv")
# m22 <- mo22 %>%
#   mutate(EARN_MDN_HI_2YR = as.numeric(EARN_MDN_HI_2YR)) %>% 
#   filter(!is.na(EARN_MDN_HI_2YR)) 






# merging data ------------------------------------------------------------

all_data <- c(
  m,
  fd1,
  fd2,
  fd3,
  fd4,
  fd5,
  fd6
  # m97,m98,m99,m00,m01,m02,m03,m04,m05,m06,m07,m08,m09,m10,m11,m12,m13,m14,
  # m15,m16,m17,m18,m19,m20,m21,m22
)

# merge datasets 

fd1 <- as.data.frame(fd1)
fd2 <- as.data.frame(fd2)
fd3 <- as.data.frame(fd3)
fd4 <- as.data.frame(fd4)
fd5 <- as.data.frame(fd5)
fd6 <- as.data.frame(fd6)
m   <- as.data.frame(m)


da <- rbind(fd1,fd2)
db <- rbind(da,fd3)
dc <- rbind(db,fd4)
dd <- rbind(dc,fd5)
de <- rbind(dd,fd6)
dg <- rbind(de,m)

combined_all_data <- dg

dartmouth_picked_data <- combined_all_data %>% 
  filter(INSTNM %in% dart_picks)
# View(dartmouth_picked_data)
# nrow(dartmouth_picked_data) #928 data pts

picked_dartmouth_data <- combined_all_data %>% 
  filter(INSTNM %in% picked_dartmouth)
# View(picked_dartmouth_data)
# nrow(picked_dartmouth_data) #1226

dartmouth_association_data <- combined_all_data %>% 
  filter(INSTNM %in% dartmouth_association) #1656 data points




# Ranking variables based on how many NAs ---------------------------------

# dartmouth picked:
# changing everything that is PS to NA
dartmouth_picked_data <- dartmouth_picked_data %>%
  mutate(across(everything(), ~ ifelse(. == "PS", NA, .)))


ranked_columns <- dartmouth_picked_data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "na_count") %>%
  arrange(na_count) %>%
  mutate(rank = row_number())
View(ranked_columns)

#____________________
#picked dartmouth
picked_dartmouth_data <- picked_dartmouth_data %>%
  mutate(across(everything(), ~ ifelse(. == "PS", NA, .)))

ranked_columns_picked_D <- picked_dartmouth_data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "na_count") %>%
  arrange(na_count) %>%
  mutate(rank = row_number())

# View the ranked columns
View(ranked_columns_picked_D)

#dartmouth assocation 
ranked_columns_dart_association <- dartmouth_association_data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "na_count") %>%
  arrange(na_count) %>%
  mutate(rank = row_number())
View(ranked_columns_dart_association)
#opeid id is institution id

#good variables:
#EARN_MDN_HI_2YR
#EARN_COUNT_WNE_HI_2YR

#EARN_CNTOVER150_HI_2Y
#EARN_COUNT_NWNE_HI_2YR

# repeating sorting process for uconn -------------------------------------
uconn_picked <- c(
  "University of Washington",
  "University of California-Los Angeles",
  "University of California-San Diego",
  "University of California-Berkeley",
  "Purdue University, West Lafayette", #seems to not have any
  "Ohio State University",
  "University of Illinois Urbana-Champaign",
  "The University of Texas at Austin",
  "University of California-Santa Barbara",
  "University of California-Davis",
  "University of Michigan-Ann Arbor",
  "University of Florida",
  "University of Wisconsin-Madison",
  "University of California-Irvine",
  "University of Georgia",
  "University of North Carolina at Chapel Hill",
  "College of William and Mary",
  "Pennsylvania State University-Main Campus",
  "University of Virginia",
  "University of Maryland-College Park",
  "Georgia Institute of Technology"
)

picked_uconn <- c(
  "North Dakota State University",
  "Arizona State University",
  "Quinnipiac University",
  "Florida State University",
  "Albany College of Pharmacy and Health Sciences",
  "SUNY at Albany",
  "Texas Tech University",
  "University of Massachusetts-Amherst",
  "West Virginia University",
  "University of New Hampshire",
  "University of Rhode Island",
  "University of Illinois Chicago",
  "Quinnipiac University",
  "University of Delaware",
  "Bentley University",
  "George Mason University",
  "The University of Texas at Arlington"
)


uconn_association<- c(picked_uconn, uconn_picked)
uconn_association <- unique(uconn_association)
print(uconn_association)

picked_uconn_data <- combined_all_data %>% 
  filter(INSTNM %in% picked_uconn)
# View(picked_uconn_data)
nrow(picked_uconn_data) 

uconn_picked_data <- combined_all_data %>% 
  filter(INSTNM %in% uconn_picked)
# View(uconn_picked_data)
nrow(uconn_picked_data) #2690


uconn_association_data <- combined_all_data %>% 
  filter(INSTNM %in% uconn_association)
#View(uconn_association_data)
nrow(uconn_association_data) #4318

#heleper function for finding unique
# unique_institutions <- unique(picked_uconn_data$INSTNM)
# unique_institutions


#function for finding proper university names
names_institutions <- combined_all_data %>%
  filter(grepl("Columbia", INSTNM)) %>%
  select(INSTNM)

names_institutions




# sorting Stanford --------------------------------------------------------
#staford picked: Johns Hopkins, Columbia, MIT, Dartmouth, Cornell, UPenn, Brown, 
# Harvard, Yale, Princeton

#picked Stanford
#NYU, Harvey Mudd College, U of Notre Dame, Washington U ST Louis, Canetgie Melon, 
#> Vanderbilt, Georgia inst of tech, Brandeis, UPenn, Yale, Cornell, USC, northwestern,
#> Art Center College of Design, Coooper Union for advancement of science and Art,
#> MIT, U of rochester, Wellesley College, Johns HOpkins, Dartmouth, Rice U, 
#> U of Michigan Ann arbor, barnard, U chicago, Harvard, Cal tech
#> 


stanford_picked <- c(
  "Johns Hopkins University", 
  "Columbia College", 
  "Massachusetts Institute of Technology", 
  "Dartmouth College", 
  "Cornell University", 
  "University of Pennsylvania", 
  "Brown University", 
  "Harvard University", 
  "Yale University", 
  "Princeton University",
  "Stanford University"
)

picked_stanford <- c(
  "Stanford University", 
  "New York University", 
  "Harvey Mudd College", 
  "University of Notre Dame", 
  "Washington University in St Louis", 
  "Carnegie Mellon University", 
  "Vanderbilt University", 
  "Georgia Institute of Technology-Main Campus", 
  "Brandeis University", 
  "University of Pennsylvania", 
  "Yale University", 
  "Cornell University", 
  "University of Southern California", 
  "Northwestern University", 
  "Art Center College of Design", 
  "The Cooper Union for the Advancement of Science and Art", 
  "Massachusetts Institute of Technology", 
  "University of Rochester", 
  "Wellesley College", 
  "Johns Hopkins University", 
  "Dartmouth College", 
  "Rice University", 
  "University of Michigan-Ann Arbor", 
  "Barnard College", 
  "University of Chicago", 
  "Harvard University", 
  "California Institute of Technology"
)

picked_stanford_data <- combined_all_data %>% 
  filter(INSTNM %in% picked_stanford)


stanford_picked_data <- combined_all_data %>% 
  filter(INSTNM %in% stanford_picked)

stanford_association<- c(picked_stanford, stanford_picked)
stanford_association <- unique(stanford_association)

stanford_association_data <- combined_all_data %>% 
  filter(INSTNM %in% stanford_association)
#View(uconn_association_data)
nrow(stanford_association_data) #2504

# dartmouth vs uconn plot -------------------------------------------------
library(ggplot2)

# Add university labels 
dartmouth_association_data$university <- "Dartmouth"
uconn_association_data$university <- "UConn"

dart_con <- rbind(dartmouth_association_data, uconn_association_data)

# Create the scatterplot
ggplot(dart_con, aes(x = university, y = EARN_MDN_HI_2YR, color = university)) +
  geom_jitter(width = 0.1) +
  labs(title = "Comparison of 2-Year Median Earnings (Dartmouth vs UConn)",
       x = "University",
       y = "2-Year Median Earnings") +
  theme_minimal()




# dartmouth vs stanford plot ----------------------------------------------
# Add university labels 
dartmouth_association_data$university <- "Dartmouth"
stanford_association_data$university <- "Stanford"

dart_stanford <- rbind(dartmouth_association_data, stanford_association_data)

# Create the scatterplot
ggplot(dart_stanford, aes(x = university, y = EARN_MDN_HI_2YR, color = university)) +
  geom_jitter(width = 0.1) +
  labs(title = "Comparison of 2-Year Median Earnings (Dartmouth vs Stanford)",
       x = "University",
       y = "2-Year Median Earnings") +
  theme_minimal()

# quantiles
quantile(dartmouth_association_data$EARN_MDN_HI_2YR, probs = c(0.25, 0.5, 0.75))
quantile(stanford_association_data$EARN_MDN_HI_2YR, probs = c(0.25, 0.5, 0.75))




x <- seq(from(min(dartmouth_association_data$EARN_MDN_HI_2YR),
              to = max(dartmouth_association_data$EARN_MDN_HI_2YR)),by = 1)
data <- dartmouth_association_data$EARN_MDN_HI_2YR
cov(dartmouth_association_data$EARN_MDN_HI_2YR, stanford_association_data$EARN_NE_MDN_3YR)

# attempt using something -------------------------------------------------

# Compute variance for all numeric columns in the dataset

# Select only numeric columns
numeric_data <- dartmouth_picked_data %>%
  select(where(is.numeric))

# Handle missing values (e.g., by removing rows with NA or imputing values)
cleaned_data <- na.omit(numeric_data)

# Compute variance-covariance matrix
variance_matrix <- var(cleaned_data)

# Check if there are any columns or rows with NA in the matrix
# Filter the matrix to remove rows and columns with NA variances
valid_indices <- complete.cases(variance_matrix)  # Identify complete rows/columns
filtered_matrix <- variance_matrix[valid_indices, valid_indices]

# 7 dartmouth data pts to look at:
# major: CiPdesc
# EARN_MDN_HI_2YR median earnings afater 2 yrs
# EARN_IN_STATE_5_YR




#ranking the dataframe based on how many variables are na
colnames(dartmouth_picked_data) %>% 
  filter(dartmouth_picked_data !=is.na)

# Convert data to a matrix for consistency
XA <- as.matrix(dartmouth_picked_data)

# Define variable names for labeling
variable_names <- colnames(XA)

# Plotting: Probability and Density
par(mfrow = c(2, 4))  # Set layout
# Add your plotting code here

for (i in 1:ncol(XA)) {
  column_data <- XA[, i]
  
  # Check if column is numeric (density) or categorical (bar plot)
  if (is.numeric(column_data)) {
    # Density plot for continuous variables
    plot(density(column_data, na.rm = TRUE),
         main = paste("Density:", variable_names[i]),
         xlab = variable_names[i],
         ylab = "Density")
    rug(column_data)
  } else {
    # Bar plot for categorical variables
    unique_vals <- sort(unique(column_data))
    proportions <- sapply(unique_vals, function(val) mean(column_data == val, na.rm = TRUE))
    barplot(proportions,
            main = paste("Probability:", variable_names[i]),
            xlab = variable_names[i],
            ylab = "Probability",
            names.arg = unique_vals)
  }
}


XA=as.matrix(dartmouth_picked_data)
n=nrow(XA)
#namit=c("High school curriculum","SAT score","College interview","Out of school activity","Sport & research programs","Essay","Letters of recommendation")
#nam=c("HSC","SAT","CI","OSA","SR","ES","LR") #short names
m=length(XA)

if(job==1) #plot pdf/pmf
{	
  m=ncol(XA)
  par(mfrow=c(2,4))
  for(i in 1:m)
  {
    xc=XA[,i]
    if(i==1 || i>4)
    {
      pr=rep(0,5)
      for(j in 1:5) pr[j]=mean(xc==j)
      plot(1:5,pr,type="h",lwd=5,ylab="Probability",xlab=nam[i],main=namit[i])		
    }
    if(i==2)
    {
      plot(density(xc),type="l",ylab="Density",xlab=nam[i],main=namit[i])
      rug(xc)				
    }
    if(i==3 || i==4)
    {
      pr=rep(0,10)
      for(j in 1:10) pr[j]=mean(xc==j)
      plot(1:10,pr,type="h",lwd=5,ylab="Probability",xlab=nam[i],main=namit[i])		
    }	
  }
}


# # View the resulting square matrix
# W = filtered_matrix
# eigenW = eigen(W,sym = T)
# p.max = eigenW$vectors[,1]
# lambda.max = eigenW$values[1]
# X = filtered_matrix
# proj1=(X-rep(1,nrow(X)))%*%t(colMeans(X))%*%p.max

#Trying again
#going t use 7 score critearia 





# Dartmouth Like Schools from data --------------------------------------------------

dart_choice <- m[m$INSTNM %in% dart_picks,]
median_earnings_2_yrs_post_grad <- dart_choice %>% 
  select(INSTNM, EARN_MDN_HI_2YR)
nrow(median_earnings_2_yrs_post_grad) #2165 data points

# converting values to numeric
median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR <- as.numeric(median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR)

# Filter out rows where EARN_MDN_HI_2YR is NA
median_earnings_no_na <- median_earnings_2_yrs_post_grad %>%
  filter(!is.na(EARN_MDN_HI_2YR))
# after filtering : nrow(median_earnings_no_na) #464



# schools that selected dartmouth:
# Johns Hopkins, Rice U, Bucknell U, Wesleyan University, Bowdoin C, Cooper Union
# For the Advancement of Science and Art, Colgate U, Tufts U, College of William 
#and Mary, Harvey Mudd College, University of Notre Dame, Wake Forest University,
# Alabama A and M university, Georgetown University, barnard College, yale, university of pennsylavia
# vanderbilt, stanford, harvard, brown



picks_dart <- m[m$INSTNM %in% picked_dartmouth,]
picks_median_earnings_2_yrs_post_grad <- picks_dart %>% 
  select(INSTNM, EARN_MDN_HI_2YR)
nrow(picks_median_earnings_2_yrs_post_grad) #2165 data points

# converting values to numeric
picks_median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR <- 
  as.numeric(picks_median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR)

# Filter out rows where EARN_MDN_HI_2YR is NA
picks_median_earnings_no_na <- picks_median_earnings_2_yrs_post_grad %>%
  filter(!is.na(EARN_MDN_HI_2YR))
# after filtering : nrow(median_earnings_no_na) #613



# Using fd1 ---------------------------------------------------------------
dart_choice_fd1 <- fd1[fd1$INSTNM %in% dart_picks,]
fd1median_earnings_2_yrs_post_grad <- dart_choice %>% 
  select(INSTNM, EARN_MDN_HI_2YR, CIPDESC)

# earn_medn_hi_2yr is the 2 years after 2 years being hired
#nrow(fd1median_earnings_2_yrs_post_grad) #2165 data points

# converting values to numeric
fd1median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR <- as.numeric(fd1median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR)

# Filter out rows where EARN_MDN_HI_2YR is NA
fd1median_earnings_no_na <- fd1median_earnings_2_yrs_post_grad %>%
  filter(!is.na(EARN_MDN_HI_2YR))

nrow(fd1median_earnings_no_na) #464
#-------------------------------------------------------------

fd1picks_dart <- fd1[fd1$INSTNM %in% picked_dartmouth,]
fd1picks_median_earnings_2_yrs_post_grad <- fd1picks_dart %>% 
  select(INSTNM, EARN_MDN_HI_2YR)
nrow(fd1picks_median_earnings_2_yrs_post_grad) #2165 data points

# converting values to numeric
fd1picks_median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR <- 
  as.numeric(fd1picks_median_earnings_2_yrs_post_grad$EARN_MDN_HI_2YR)

# Filter out rows where EARN_MDN_HI_2YR is NA
fd1picks_median_earnings_no_na <- fd1picks_median_earnings_2_yrs_post_grad %>%
  filter(!is.na(EARN_MDN_HI_2YR))
# after filtering : nrow(fd1picksmedian_earnings_no_na) 613



# checks on the dataset ---------------------------------------------------
# # ensures it contains all the names I was hoping for 
# median_earnings_2_yrs_post_grad %>% 
#   distinct(INSTNM)

# Dartmouth Data  ---------------------------------------------------------
# dartmouth <- m[m$INSTNM == "Dartmouth College",]
# dartmouth_undergrad <- dartmouth[dartmouth$CREDDESC== "Bachelor's Degree",]
# #51 dartmouth undergrad data pts to work with
# #most earnings info is surpressed
# 
# dartmouth_grad <- dartmouth[dartmouth$CREDDESC == "Master's Degree",]
