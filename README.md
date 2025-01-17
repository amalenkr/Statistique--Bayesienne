# Prediction d'accident d'avion

Dans le cadre de notre projet de statistique bayésienne, nous avons travaillé sur la prédiction du nombre moyen d’accidents d’avion à partir de la base de donnée "plane_accident". 

Chaque ligne de cette base de données correspond à un accident d’avion ayant eu lieu entre le 01/01/1972 et le 31/12/1975 et contient 3 colonnes : 
-	La date exacte de l’accident
-	Le nombre de victimes
-	Le jour de l’accident (allant de 1 -> 01/01/1972 à 1460=4*365 -> 31/12/1975).

Nous avons comparé une approche visant à estimer le nombre d’accidents moyen par semaine et par mois grâce à une non-informative prior et une informative prior.

Nous avons adopté l’approche suivante :
1)	Justifier la pertinence d’utiliser une loi de poisson pour modéliser l’apparition d’évènements dans le temps (des accidents d’avion dans notre cas)
2)	Calculer la loi a posteriori de théta lorsque X suit une loi de poisson de moyenne théta et que la loi a priori de théta est une gamma (informative)
3)	Calculer la loi a posteriori de théta lorsque X suit une loi de poisson de moyenne théta et que la loi a priori de théta est une uniforme sur [0 ; max(nombre d’accident par semaine / mois)] (non-informative)
4)	Comparer les 2 distributions empiriques par semaines & par mois
5)	Enfin réaliser un Bootstrap de nos observations pour comparer nos estimateurs bayésiens


