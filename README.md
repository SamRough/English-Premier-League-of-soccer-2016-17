# English Premier League of soccer(2016-17)

**Download the file EPL_1617.csv, which contains details on each game played in the 2016/17 season of the English Premier League of soccer, including the date, home team, away team and number of goals scored by each team.**

1. Construct an R function named `season.summary` that produces, given the outcomes to all the games of a soccer season, a table as the one shown in Figure 1 below. The returned table should be an integer matrix with appropriate column names and row names corresponding to the team names for improved readability.<br/>Note that the number of points accumulated by each team are obtained using the following very simple rule:
	* 3 points for each win, 
	* 1 point for each tie,
	* no points for losses.<br/>*Note also that teams are ordered*
	* first according to total points,
	* then, in case of a tie, in terms of goal difference,
	* finally, if two teams are still tied, in terms of the number of wins.

2. Construct a function named `dated.summary` that returns a table like the one produced in Assignment 1, but only accounting for games that took place between two dates, first.date and second.date, to be provided by the user as arguments. These dates can be expected to be character strings in the format ’yyyy-mm-dd’. The function should also use the game.results for the whole season (you will be passing on the soccer data) as an argument.

3. Construct an R function named `team.progression` that produces, given the game.results, a table showing how the points total of each team progressed over the season after each game,Specifically, the output of the function should be a matrix with appropriate team names as row names for improved readability.

4. graph that can be used to assess Chelsea’s domination over the season. The specifics of what the graph should look like are up to you, but you should keep in mind that a lot of this is about how to display the information clearly and in a way that is insightful. One approach would be to compare Chelsea’s performance to that of a few carefully selected teams. Write a few comments to explain what conclusion you draw from your graph. Finally, do the same to assess Sunderland’s performance over the season.

5. A **Poisson regression model** for English Premiership Soccer. To model the scores of soccer games from the English Premier League, assume that in every league game, 
	* the number of goals X scored by team i when playing team j at home satisfies:<br/>**X ∼ Poisson(exp(μ + ∆ + αi + βj))**
	* and that the number of goals scored by team j on the road playing team i satisfies:<br/>**Y ∼ Poisson(exp(μ + αj + βi))**
		* μ : overall league scoring parameter,
		* ∆ : scoring effect for playing at home (same for all teams in the league),
		* αi : offensive effect for team i (a large positive value indicates a strong offense),
		* βj : defensive effect for team j (a large negative value indicates a strong defense).
	
	* Create an R function `new.season` that **simulates**, based on the above Poisson model, a new season from a calendar of games, a data frame containing the dates of all the games and which **teams play home and away for each game, and the parameters mu, Delta, alpha and beta of the model**.

	* Using **the parametric bootstrap approach**, create an R function that returns **simultaneous confidence intervals** for all parameters of the model.

6. Create an R function that, given the parameters mu, Delta, alpha and beta of the model and a number of seasons N.seasons to be simulated, will
	* simulate N.seasons different seasons,
	* for each simulated season, assign a rank to each team (calculated using the same rule as
in Assignment 2),
	* return a table showing the number of times each team was awarded each rank and the average rank of each team.
<br/>Then, use your function to simulate 100 seasons on the EPL and provide the appropriate summary of your simulation.

7. At the end of a season, the **bottom three teams** in the league rankings are relegated from the EPL to a lower division of English professional soccer for the next season. Simultaneously, three new teams from lower divisions are promoted to the EPL.how likely was it for each of these teams to be relegated? Also, among the teams that were not relegated, which team was the most likely to have been relegated?

8. Finally, how likely was **Chelsea** to end up **League Champion**? Which team was most likely to finish first in the standings?there is a clear discrepancy between the top seven teams and the bottom 13 teams. The see this, and again using your previous output, calculate the approximate probability that each team finishes in the top 7.
