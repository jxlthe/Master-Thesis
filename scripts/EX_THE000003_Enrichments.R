# pyhper("success in sample - 1",
#        "success in background",
#        "failure in background",
#        "sample size",
#        lower.tail = FALSE)


set <- c()
p <- c()

#All

{
#14771  Conserved
#3768   Conserved, No DEG in C4, At least One of Each C3 DEG
#183    Conserved, All C4 Expanded
#36    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG")
  
p <- append(p, phyper(36 -1,
       3768,
       14771 - 3768,
       183,
       lower.tail = FALSE))

#14771  Conserved
#1208   Conserved, No DEG in C4, At least One of Each C3 DEG, Only Upregulated
#183    Conserved, All C4 Expanded
#14    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Upregulated

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Upregulated")

p <- append(p, phyper(14 -1,
       1208,
       14771 - 1208,
       183,
       lower.tail = FALSE))

#14771  Conserved
#915   Conserved, No DEG in C4, At least One of Each C3 DEG, Only Downregulated
#183    Conserved, All C4 Expanded
#10    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Downregulated

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Downregulated")

p <- append(p, phyper(10 -1,
       915,
       14771 - 915,
       183,
       lower.tail = FALSE))

}

#14771  Conserved
#180 Conserved At least one DEG in Zm, No C3 DE
#3786 Conserved Zm Expanded
#52 Conerved Zm Expanded, At least one DEG in Zm, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm, No C3 DE")

p <- append(p, phyper(52 - 1,
       180,
       14771 - 180,
       3786,
       lower.tail = FALSE))

#14771  Conserved
#90 Conserved At least one DEG in Zm Only Down, No C3 DE
#3786 Conserved Zm Expanded
#22 Conerved Zm Expanded, At least one DEG in Zm Only Down, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm Only Down, No C3 DE")

p <- append(p, phyper(22 - 1,
       90,
       14771 - 90,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#90 Conserved At least one DEG in Zm Only Up, No C3 DE
#3786 Conserved Zm Expanded
#30 Conerved Zm Expanded, At least one DEG in Zm Only Up, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm Only Up, No C3 DE")

p <- append(p, phyper(30 - 1,
       90,
       14771 - 90,
       3786,
       lower.tail = FALSE))


#Constitutive Hypothesis

#Zm
{

#14771  Conserved
#3937 Conserved No DEG in Zm, At least one of Each C3 DE
#3786 Conserved Zm Expanded
#974 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 DE

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 DE")
  
  p <- append(p, phyper(974 - 1,
       3937,
       14771 - 3937,
       3786,
       lower.tail = FALSE))

#14771  Conserved
#964 Conserved No DEG in Zm, At least one of Each C3 Down
#3786 Conserved Zm Expanded
#236 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Down

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Down")
  
  p <- append(p, phyper(236 - 1,
       964,
       14771 - 964,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#1258 Conserved No DEG in Zm, At least one of Each C3 Up
#3786 Conserved Zm Expanded
#345 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Up

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Up")
  
  p <- append(p, phyper(345 - 1,
       1258,
       14771 - 1258,
       3786,
       lower.tail = FALSE))


#*
#14771  Conserved
#1680 Conserved No DEG in Zm, All C3 DE
#3786 Conserved Zm Expanded
#513 Conserved Zm Expanded, No DEG in Zm, All C3 DE

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE")
  
  p <- append(p, phyper(513 - 1,
       1680,
       14771 - 1680,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#382 Conserved No DEG in Zm, All C3 DE Only Down
#3786 Conserved Zm Expanded
#105 Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Down

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Down")
  
  p <- append(p, phyper(105 - 1,
       382,
       14771 - 382,
       3786,
       lower.tail = FALSE))

#*
#14771  Conserved
#610 Conserved No DEG in Zm, All C3 DE Only Up
#3786 Conserved Zm Expanded
#204 Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Up

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Up")
  
  p <- append(p, phyper(204 - 1,
       610,
       14771 - 610,
       3786,
       lower.tail = FALSE))


}

#Ms
{
  
  #14771  Conserved
  #4812 Conserved No DEG in Ms, At least one of Each C3 DE
  #452 Conserved Ms Expanded
  #122 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 DE
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 DE")
  
  p <- append(p, phyper(122 - 1,
         4812,
         14771 - 4812,
         452,
         lower.tail = FALSE))
  
  #14771  Conserved
  #1192 Conserved No DEG in Ms, At least one of Each C3 Down
  #452 Conserved Ms Expanded
  #34 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Down
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Down")
  
  p <- append(p, phyper(34 - 1,
         1192,
         14771 - 1192,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #1560 Conserved No DEG in Ms, At least one of Each C3 Up
  #452 Conserved Ms Expanded
  #40 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Up
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Up")
  
  p <- append(p, phyper(40 - 1,
         1560,
         14771 - 1560,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #2086 Conserved No DEG in Ms, All C3 DE
  #452 Conserved Ms Expanded
  #62 Conserved Ms Expanded, No DEG in Ms, All C3 DE
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE")
  
  p <- append(p, phyper(62 - 1,
         2086,
         14771 - 2086,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #473 Conserved No DEG in Ms, All C3 DE Only Down
  #452 Conserved Ms Expanded
  #16 Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Down
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Down")
  
  p <- append(p, phyper(16 - 1,
         473,
         14771 - 473,
         452,
         lower.tail = FALSE))
  
  #14771  Conserved
  #787 Conserved No DEG in Ms, All C3 DE Only Up
  #452 Conserved Ms Expanded
  #24 Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Up
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Up")
  
  p <- append(p, phyper(24 - 1,
         787,
         14771 - 787,
         452,
         lower.tail = FALSE))
  
  
}

#Sb
{
  
  #14771  Conserved
  #4913 Conserved No DEG in Sb, At least one of Each C3 DE
  #981 Conserved Sb Expanded
  #292 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 DE
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 DE")
  
  p <- append(p, phyper(292 - 1,
         4913,
         14771 - 4913,
         981,
         lower.tail = FALSE))
  
  #14771  Conserved
  #1233 Conserved No DEG in Sb, At least one of Each C3 Down
  #981 Conserved Sb Expanded
  #80 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Down
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Down")
  
  p <- append(p, phyper(80 - 1,
         1233,
         14771 - 1233,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #1577 Conserved No DEG in Sb, At least one of Each C3 Up
  #981 Conserved Sb Expanded
  #104 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Up
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Up")
  
  p <- append(p, phyper(104 - 1,
         1577,
         14771 - 1577,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #2146 Conserved No DEG in Sb, All C3 DE
  #981 Conserved Sb Expanded
  #139 Conserved Sb Expanded, No DEG in Sb, All C3 DE
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE")
  
  p <- append(p, phyper(139 - 1,
         2146,
         14771 - 2146,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #489 Conserved No DEG in Sb, All C3 DE Only Down
  #981 Conserved Sb Expanded
  #31 Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Down
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Down")
  
  p <- append(p, phyper(31 - 1,
         489,
         14771 - 489,
         981,
         lower.tail = FALSE))
  
  #14771  Conserved
  #814 Conserved No DEG in Sb, All C3 DE Only Up
  #981 Conserved Sb Expanded
  #55 Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Up
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Up")
  
  p <- append(p, phyper(55 - 1,
         814,
         14771 - 814,
         981,
         lower.tail = FALSE))
  
  
}


results <- data.frame(set, p)
results$padj <- p.adjust(results$p, method="bonferroni")
