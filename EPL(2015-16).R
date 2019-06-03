## 项目主题
基于16/17赛季英超380场比赛的主客场球队的进球数据，做一些相应数据分析。

### 英格兰超级联赛背景知识
英超共有 20 支参赛球队，赛季于每年8月至次年5月进行，每赛季英超共有380场比赛，每支球队主场和客场比赛各有19场。每场胜方可得3分，平局各得1分，负方得0分，按各队于联赛所得的积分排列。完成所有赛事后总积分最高的队伍可以夺得联赛冠军，而总积分最低的3队球队会降级至英冠联赛。

### 话不多说，正式开始呐。

#### 初探数据
- 发现这个是一个380*5的表格，很干净的表格不需要额外的清洗整理数据。
- 5个列分别表示“日期”、“HT-主场球队”、“AT-客场球队”、“HG-主场球队进球数”、“AG-客场球队进球数”
- 380行表示380场比赛，以第一行为例，2016-08-14 16:00这个时间，主队阿森纳，客队利物浦，主队进球数3，客队进球数4。所以总结起来为主队阿森纳以比分3：4败于客队利物浦。
```{r}
# 设置数据读取路径
setwd("~/Desktop/University/UM/STAT4600/Assignment")
# 读取 "CSV" 文件
EPl.1617 <- read.csv("EPL_1617.csv", header=T, stringsAsFactors = FALSE)
# 表格维度
dim(EPl.1617)
# 前几个行的数据
head(EPl.1617)
```


#### 由于380场比赛数据在手，那么什么是我们最想要知道的事情呢。自然而然16/17赛季谁是冠军，谁被降级，这是我们最关心的事情。
- 根据英超赛制，每场胜方可得3分，平局各得1分，负方得0分
- 按照积分高低排名，如果积分一样，按照净胜球数排名，如果净胜球数也一样，按照胜场排序
- 我们想要输出一个汇总表格，显示20支队伍名称，胜场数，平场数，输场数，总积分，总进球数，总丢球数，净胜球数

#### 创建一个叫`season.summary`的function，当我们输入比分数据，来自动生成一个赛季汇总信息表格。
#### 代码如下：
```{r}
# removing the first column "date".
EPl.1617 <- EPl.1617[1:380,2:5]

# "a" is variable of function "season.summary"
season.summary <- function(a){
  wins.team <- c()
  losses.team <- c()
  ties.team <- c()
  # using "for" loop to calculate frequences of wins,losses and ties for every teams.
  for (i in 1:nrow(a)){
    if(a[i,3] > a[i,4]){
      wins.team <- c(wins.team,as.character(a[i,1]))
      losses.team <- c(losses.team,as.character(EPl.1617[i,2]))
    }
    if(a[i,3] < a[i,4]){
      wins.team <- c(wins.team,as.character(EPl.1617[i,2]))
      losses.team <- c(losses.team,as.character(a[i,1]))
    }
    if(a[i,3] == a[i,4]){
      ties.team <- c(ties.team,as.character(a[i,1]),as.character(a[i,2]))
    }
  }
  
  # number of wins <- n.wins
  # number of ties <- n.ties
  # number of losser <- n.losses
  # pts.season is a matrix of points for every team.
  # s.season is a part of final output matrix but without order.
  n.wins <- as.matrix(table(wins.team))
  n.ties <- as.matrix(table(ties.team))
  n.losses <- as.matrix(table(losses.team))
  pts.season <- 3*n.wins+1*n.ties
  s.season <- cbind(n.wins,n.ties,n.losses,pts.season)
  
  library(plyr)
  #goal scores at home
  HG.season <- count(EPl.1617,"HT","HG")
  #goal scores away home
  AG.season <- count(EPl.1617,"AT","AG")
  #total goal scores for each team
  TG.season <- cbind(HG.season,AG.season)
  TG.season <- cbind.data.frame(TG.season[,1],TG.season[,2]+TG.season[,4])
  
  
  #loss scores at home
  HL.season <- count(EPl.1617,"HT","AG")
  #loss scores away home
  AL.season <- count(EPl.1617,"AT","HG")
  #total loss scores for each team
  TL.season <- cbind(HL.season,AL.season)
  TL.season <- cbind.data.frame(TL.season[,1],TL.season[,2]+TL.season[,4])
  
  #goal differnece
  GD.season <- cbind.data.frame(TG.season[,1],TG.season[,2]-TL.season[,2])
  
  #g.season is another part of final answer matrix but without order.
  g.season <- cbind(TG.season,TL.season,GD.season)
  g.season <- g.season[,c(1:2,4,6)]
  
  #combine s.season and g.season to make up final matrix final.season but without order.
  final.season <- cbind(s.season,g.season[,2:4])
  colnames(final.season) <- c("W","T","L","Pts","GF","GA","GD")
  final.season <- as.matrix(final.season)
  
  
  # re-order the table according to decreasing point totals.
  final.season <- final.season[order(final.season[,4],final.season[,7],final.season[,1],decreasing = TRUE),]
  
  
  return(final.season)
  
}

season.summary(a = EPl.1617)
```
####基于上面输出的赛季汇总的表格，我们可以知道：
- 切尔西，托特纳姆热刺，曼城为前三强队伍，这三只球队也是仅有的三支总进球数达到80个的球队。有趣的是冠军切尔西净胜球数却比亚军热刺少了8个，可见热刺平局太多影响了总积分。
- 赫尔城，米德尔斯堡，桑德兰是排在最末尾的三支球队，他们面临降级窘境


####当有了赛季的summary数据，我们也想要知道具体在某个时间段，这些队伍的表现如何。例如，在赛季开始前的前7周哪支队伍势头最强，或者赛季结束前8周，哪支队伍表现不尽人意。
#### 创建一个叫`date.summary`的Function，当我们输入具体的’yyyy-mm-dd’时间段，来自动生成一个具体时间段内汇总信息表格。
#### 代码如下：
```{r}

```
####上面的第一个season summary调用的时候居然跑不出来，痛苦啊。先用一下Alex教授的season.summary。
```{r}
season.summary <- function(scores, alphabetical=T){
  teams <- sort(unique(c(scores$HT,scores$AT)))     # this works even if a team has only home/away games
  n.teams <- length(teams)
  n.games <- length(scores$HT)
  ##
  # Calculating W-T-L
  #
  WTL <- matrix(0,n.teams,7,dimnames=list(teams,c('W','T','L','GF','GA','GD','Pts')))
  for(j in 1:n.games){
    home <- scores$HT[j]
    away <- scores$AT[j]
    WTL[home,'GF'] <- WTL[home,'GF'] + scores$HG[j]
    WTL[home,'GA'] <- WTL[home,'GA'] + scores$AG[j]
    WTL[away,'GF'] <- WTL[away,'GF'] + scores$AG[j]
    WTL[away,'GA'] <- WTL[away,'GA'] + scores$HG[j]
    if(scores$HG[j] > scores$AG[j]) outcome <- 'H'          # home win in this case 
    if(scores$HG[j] < scores$AG[j]) outcome <- 'A'          # away win
    if(scores$HG[j] == scores$AG[j]) outcome <- 'T'         # Tie
    switch(outcome,
           'A' = {WTL[home,'L'] <- WTL[home,'L']+1 ; WTL[away,'W'] <- WTL[away,'W']+1},
           'T' = {WTL[home,'T'] <- WTL[home,'T']+1 ; WTL[away,'T'] <- WTL[away,'T']+1},
           'H' = {WTL[home,'W'] <- WTL[home,'W']+1 ; WTL[away,'L'] <- WTL[away,'L']+1}
    )
  }
  #
  # Calculating goal difference and points
  #
  WTL[,'GD'] <- WTL[,'GF']-WTL[,'GA']
  WTL[,'Pts'] <- 3*WTL[,'W']+WTL[,'T']
  #
  # In alphabetical order
  #
  if(alphabetical==T) return(WTL)
  #
  # Otherwise, by points total, GD and Wins
  #
  return(WTL[order(WTL[,'Pts'],WTL[,'GD'],WTL[,'W'],decreasing=T),])
}

```

####代码如下
```{r}
date.summary <- function(first.date, second.date, game.results){
  first.date <- as.Date(first.date)
  second.date <- as.Date(second.date)
  game.results.Date <- as.Date(game.results$Date)
  
  if(first.date - second.date >= 0){
    return(paste("The two dates aren???t ordered properly."))
  }
  if(max(game.results.Date) <= first.date | min(game.results.Date) >= second.date){
    return(paste("no games took place between the given dates"))
  } 
  if(ncol(game.results[game.results$Date >= first.date & game.results$Date <= second.date,]) > 0){
    game.results <- game.results[game.results$Date >= first.date & game.results$Date <= second.date,]
    return(season.summary(game.results))
  }
}
```
####看看"2016-08-13","2016-09-30"和"2017-04-03","2017-05-21"两个时间段内的表现如何。
```{r}
setwd("~/Desktop/University/UM/STAT4600/Assignment")
game.results <- read.csv("EPL_1617.csv", header = T, stringsAsFactors = FALSE)
# To find the first seven weeks of the season. "2016-08-13" to "2016-09-30"
# In the first seven weeks "Manchester City" is the best. "Sunderland" is the worst.
date.summary("2016-08-13","2016-09-30",game.results)
# To find the last seven weeks of the season. "2017-04-03" to "2017-05-21"
# In last first seven weeks "Tottenham Hotspur" is the best. "West Bromwich Albion" is the worst.
date.summary("2017-04-03","2017-05-21",game.results)

```


####到现在为止我们已经知道每支球队的总积分，排名，总进球数等信息，我们已经知道了这些结果，但是我们也想要知道这中间发生了什么，每支球队赛季中具体的表现如何，球队是个什么样的走势。

####首先定义一个function`team.progression`，这里progression表示积分累积的意思。当我们输入“比分数据”进这个function，它就会输出每支球队积分的累积过程。
```{r}
team.progression <- function(game.results){
  k <- rep(0,20)
  count.pts <- rep(0,20)
  teams <- game.results[order(game.results[,1]),]
  teams.name <- sort(unique(c(game.results$HT,game.results$AT)))
  result.matrix <- matrix(0,20,38,byrow = T)
  colnames(result.matrix) <- as.character(c(1:38)) 
  rownames(result.matrix) <- teams.name
  for(i in 1:20){
    for(j in 1:380){
      if(teams.name[i] == teams[j,2]){
        k[i] <- k[i]+1
        if(teams[j,4] > teams[j,5]){
          count.pts[i] <- count.pts[i] + 3
        }
        if(teams[j,4] < teams[j,5]){
          count.pts[i] <- count.pts[i] + 0
        }
        if(teams[j,4] == teams[j,5]){
          count.pts[i] <- count.pts[i] + 1
        }
        result.matrix[i,k] <- count.pts[i]
      }
      if(teams.name[i] == teams[j,3]){
        k[i] <- k[i]+1
        if(teams[j,4] > teams[j,5]){
          count.pts[i] <- count.pts[i] + 0
        }
        if(teams[j,4] < teams[j,5]){
          count.pts[i] <- count.pts[i] + 3
        }
        if(teams[j,4] == teams[j,5]){
          count.pts[i] <- count.pts[i] + 1
        }
        result.matrix[i,k] <- count.pts[i]
      }
    }
  }
  
  return(result.matrix)
}

team.progression(game.results)

```
#### 数据不够直观，别担心，我们可以把它画出来
#### 首先看看“切尔西”表现如何
```{r}
# Chelsea's domination over the season.
for(i in 1:20){
  plot(team.progression(game.results)[i,],xlim = c(0,40),ylim = c(0,100),xlab="",ylab = "",type = "l")
  par( new = TRUE)
}
par(new = T)
plot(team.progression(game.results)[4,],xlim = c(0,40),ylim = c(0,100),xlab="game 1-38",ylab = "points",type = "l",col = "red")
title(main = "points over the season")
```
####很明显可以看出赛季冠军切尔西，在图中一枝独秀，虽然赛季初表现一般，不过大概十场以后开始发力，并一直保持优势到38场比赛结束。
####再看看最后一名“桑德兰"同学的表现如何。
```{r}
# Sunderland's performance over the season.
for(i in 1:20){
  plot(team.progression(game.results)[i,],xlim = c(0,40),ylim = c(0,100),xlab="",ylab = "",type = "l")
  par( new = TRUE)
}
par(new = T)
plot(team.progression(game.results)[15,],xlim = c(0,40),ylim = c(0,100),xlab="game 1-38",ylab = "points",type = "l",col = "green")
title(main = "points over the season")
```


### 如果想基于现有的数据，我们知道了这些排名，通过排名可以反映实例的强弱，但是我们也想量化每个球队的进攻能力，防守能力，主场优势等等。首先我们明确，主场和客场进球数是属于泊松分布的。基于Lee的论文，我们得到这么一个模型。这里每个参数都有相对应的解释。
#### mu表示          delta表示        alpha表示       beta表示

#### 第一步，找到各个球队对应的参数，这里我们用MLE来估计参数。因为因变量进球数是离散型的变量，所有我们不用lm，改用glm做简单的线性回归，得到估计参数
```{r}
soccer.estimates <- function(scores){
  #
  # Building the X matrix for GLM fit
  #
  home.teams <- as.character(scores[,'HT'])
  away.teams <- as.character(scores[,'AT'])
  teams <- sort(unique(c(home.teams,away.teams)))   
  n.teams <- length(teams)
  n.games <- length(home.teams)
  Home <- matrix(0,n.games,n.teams-1)
  Away <- matrix(0,n.games,n.teams-1)
  for(j in 1:n.games){
    home <- which(home.teams[j]==teams)
    away <- which(away.teams[j]==teams)
    if(home <= (n.teams-1)){
      Home[j,] <- rep(0,n.teams-1)
      Home[j,home] <- 1
    }
    else{
      Home[j,] <- rep(-1,n.teams-1)
    }
    if(away <= (n.teams-1)){
      Away[j,] <- rep(0,n.teams-1)
      Away[j,away] <- 1
    }
    else{
      Away[j,] <- rep(-1,n.teams-1)
    }
  }
  home.matrix <- cbind(rep(1,n.games),rep(1,n.games),Home,Away)
  away.matrix <- cbind(rep(1,n.games),rep(0,n.games),Away,Home)
  design.matrix<-rbind(home.matrix,away.matrix)
  #
  # The response vector
  #
  home.scores <- scores[,'HG']
  away.scores <- scores[,'AG']
  response.vector <- c(home.scores,away.scores)
  #
  # Fitting the Poisson model
  #
  model.fit <- glm(response.vector~0+design.matrix,family=poisson(link='log'))
  mu <- model.fit$coef[1]
  names(mu) <- NULL
  delta <- model.fit$coef[2]
  names(delta) <- NULL
  alpha <- c(model.fit$coef[3:21],-sum(model.fit$coef[3:21]))
  names(alpha) <- teams
  beta <- c(model.fit$coef[22:40],-sum(model.fit$coef[22:40]))
  names(beta) <- teams
  return(list(mu=mu,delta=delta,alpha=alpha,beta=beta))
}

specific.par <- soccer.estimates(game.results)
specific.par
```

#### 基于上面算出来的参数，我们可以通过这些参数来预测新一轮的进球数。
```{r}
############################################################
# define a function new season.

new.season <- function(calendar,mu,delta,alpha,beta){
  
  home.teams <- as.character(calendar[,'HT'])
  away.teams <- as.character(calendar[,'AT'])
  teams <- sort(unique(c(home.teams,away.teams)))   
  par.mu <- rep(mu,nrow(calendar))
  par.delta <- rep(delta, nrow(calendar))
  par.alpha.i <- c(0,nrow(calendar))
  par.beta.i <- c(0,nrow(calendar))
  par.alpha.j <- c(0,nrow(calendar))
  par.beta.j <- c(0,nrow(calendar))
  
  for(i in 1:length(teams)){
    a <- which(home.teams == teams[i])
    par.alpha.i[a] <- alpha[teams[i],]
    par.beta.i[a] <- beta[teams[i],]
  }
  
  for(i in 1:length(teams)){
    a <- which(away.teams == teams[i])
    par.alpha.j[a] <- alpha[teams[i],]
    par.beta.j[a] <- beta[teams[i],]
  }
  # table1 gives a data.frame containing calendar of games and the parameters mu, Delta, alpha and beta of the model.
  table1 <- cbind(calendar,par.mu,par.delta,par.alpha.i,par.beta.i,par.alpha.j,par.beta.j)
  
  # table2 gives a data frame containing the calendar and scores to all the games 
  # in the same format as the soccer object we have used before.
  HG <- apply(table1[,4:9],1,function(x) rpois(n=1,lambda = exp(x[1]+x[2]+x[3]+x[6])))
  AG <- apply(table1[,4:9],1,function(x) rpois(n=1,lambda = exp(x[1]+x[5]+x[4])))
  table2 <- cbind(calendar, HG, AG)
  
  return(table2)
}
setwd("~/Desktop/University/UM/STAT4600/Assignment")
soccer <- read.csv("EPL_1617.csv", header=T, stringsAsFactors = FALSE)

estimates <- soccer.estimates(soccer)
calendar <- game.results[,1:3]
mu <- as.numeric(estimates$mu)
delta <- as.numeric(estimates$delta)
alpha <- as.data.frame(estimates$alpha)
beta <- as.data.frame(estimates$beta)

#producing a new season.
NewSeason <- new.season(calendar,mu,delta,alpha,beta)
head(NewSeason)
```

####我们预测了新一个赛季的比分，基于新一个赛季的比分，下面用 the parametric bootstrap approach来simultaneous N个赛季的parameters

```{r}
# the game.results of a season of interest
# the number N of replicates of the full season that should be used for bootstrapping
# the confidence level C that is desired

# define a par.bootstrap function to find confidence interval for all parameters.
par.bootstrap <- function(calendar,mu,delta,alpha,beta, N, C){
  teams <- sort(unique(c(as.character(EPl.1617[,'HT']),as.character(EPl.1617[,'AT'])))) 
  # for every bootstrap, get 42 pars and put it in a single cols. N times bootstrap, N cols pars.
  bootstrap.sample <- matrix(NA,nrow = 42,ncol = N)
  for(i in 1:N){
    bootstrap.sample[,i] <- unlist(soccer.estimates(new.season(calendar,mu,delta,alpha,beta)), use.names = FALSE)
  }
  CI <- t(apply(bootstrap.sample,1,function(x) quantile(x,probs=c(C/2,1-C/2))))
  name.alpha <- paste(teams, "alpha",sep = "__")
  name.beta <- paste(teams, "beta",sep = "__")
  row.names(CI) <- c('mu','Delta',name.alpha,name.beta)
  return(CI)
}
# if N is large, may take a while to get answer.
CI.bootstrap <- par.bootstrap(calendar,mu,delta,alpha,beta,N=100,C=0.1)
CI.bootstrap

```
#### 这样我们就量化了每个球队的进攻和防守强弱，主场优势等数值。

#### 然后我们看看平均下来，各个球队的排名如何
```{r}
# First defining a function that returns the ranks of each team at the end
# of a season, given the scores to all the games (using the summary function
# of Assignment 2)
#
season.ranks <- function(scores){
  summary.table <- season.summary(scores,F)			# teams are then ranked
  ranks <- order(rownames(summary.table))
  return(ranks)
}
#
# A function for calculating the average rank from a vector of rank counts
#
avgrank <- function(obs.ranks){
  sum((1:length(obs.ranks))*obs.ranks)/sum(obs.ranks)
}

# Simulating seasons and producing rank summaries
#
simulate.seasons <- function(calendar,mu,Delta,alpha,beta,N.seasons){
  teams <- row.names(alpha)
  n.teams <- length(teams)
  team.ranks <- matrix(0,n.teams,n.teams,dimnames=list(teams,1:n.teams))		# the table to be returned
  for(j in 1:N.seasons){
    new.ranks <- season.ranks(new.season(calendar,mu,Delta,alpha,beta))		# the ranks for season j
    current.ranks <- matrix(0,n.teams,n.teams)								# the matrix of rank indicators
    for(k in 1:n.teams){
      current.ranks[k,new.ranks[k]] <- 1
    }
    team.ranks <- team.ranks+current.ranks			# adds one instance of appropriate rank for each team
  }
  average.rank <- apply(team.ranks,1,avgrank)						# the average rank for each team
  return(cbind(team.ranks, AvgRank=average.rank))
}
#
rank.table <- simulate.seasons(calendar, mu, delta, alpha, beta, 100)
rank.table
```
#### 排名1-20，我们simulate了100个赛季，切尔西28次夺冠，热刺61次，利物浦4次曼城5次，曼联2次

#### 看看每支球队降级的概率：
```{r}
#
# To obtain the approximate probabilities of relegation for all the teams
#
rowSums(rank.table[,18:20])/100
```

####切尔西夺冠的概率是多少
```{r}
#
# From reading the first column!
#
rank.table[,1]/100
```

```{r}
#
# For Top 7 finishes
#
rowSums(rank.table[,1:7])/100
```
#### 阿森纳、切尔西、曼城、热刺、曼联和利物浦这几次球队几乎100%会进入前7。正好他们也是几支英超传统强队。

