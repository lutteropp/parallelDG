# parallelDG

Dies ist eine Ausarbeitung für das Abschlussprojekt des Praktikums Multicore-Programmierung
im Wintersemester 2015/16. Ziel des Projektes war es, am Beispiel des Jacobi-Verfahrens und
des Gauß-Seidel-Verfahrens parallele Lösungsmethoden partieller Differentialgleichungen zu
implementieren. Wir haben sowohl das Jacobi-Verfahren als auch das Gauß-Seidel-Verfahren
mittels OpenMP parallelisiert. Hierbei haben wir für das Gauß-Seidel-Verfahren sowohl einen
Parallelisierungsansatz mittels Wavefront als auch einen Parallelisierungsansatz mittels Rot-
Schwarz-Unterteilung gewählt. Unsere mit SIMD-Instruktionen beschleunigte Rot-Schwarz-
Implementierung des Gauß-Seidel-Verfahrens zeigte die schnellste Laufzeit bei unseren Messungen, gefolgt von der nicht vektorisierten Rot-Schwarz-Implementierung.

Unsere Ausarbeitung ist in doc/doku.pdf zu finden.
