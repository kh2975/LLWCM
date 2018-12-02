install.packages("dplyr")
library(dplyr)
find_candidates = function(file){
  
  
  mydata <- read.csv(file, header=TRUE, sep=",")

  #return(mydata)
  
  MS_lengths <- select(mydata, ms_start, ms_end,Length.of.MS.repeat)
  
  sorted1 <- mydata[with(mydata, order(-STAR_cells_in_gene, ms_distance_from_3prime, -Length.of.MS.repeat)), ]
  #sorted11 <- mydata[with(mydata, order(-STAR_cells_in_gene, -Length.of.MS.repeat, ms_distance_from_3prime)), ]
  sorted_by_distance <- mydata[with(mydata, order(ms_distance_from_3prime, -Length.of.MS.repeat,-STAR_cells_in_gene)), ]
  #sorted22 <- mydata[with(mydata, order(ms_distance_from_3prime,-STAR_cells_in_gene,-Length.of.MS.repeat)), ]
  sorted3 <- mydata[with(mydata, order(-Length.of.MS.repeat, -STAR_cells_in_gene, ms_distance_from_3prime)), ]
  #sorted33 <- mydata[with(mydata, order(-Length.of.MS.repeat, ms_distance_from_3prime, -STAR_cells_in_gene)), ]
  
  sorted1<- select(sorted1,specific_gene_ncbi_id, ms_start, ms_end,STAR_cells_in_gene, ms_distance_from_3prime,Length.of.MS.repeat)
  #s2<- select(sorted2,specific_gene_ncbi_id, ms_start, ms_end,STAR_cells_in_gene, ms_distance_from_3prime,Length.of.MS.repeat)
  sorted3<- select(sorted3,specific_gene_ncbi_id, ms_start, ms_end,STAR_cells_in_gene, ms_distance_from_3prime,Length.of.MS.repeat)
  
  distance_no_mono <- with(sorted_by_distance, which(repeat_type=="mono", arr.ind=TRUE))
  sorted_distance_no_mono <- sorted_by_distance[-distance_no_mono, ]
  sorted_distance_with_mono <- sorted_by_distance[distance_no_mono, ]
  
  highlycovered_nomono <- with(sorted_distance_no_mono, which(STAR_cells_in_gene<3454, arr.ind=TRUE))
  sorted_distance_nomono_highlycovered <- sorted_distance_no_mono[-highlycovered_nomono, ]
  
  highlycovered_mono <- with(sorted_distance_with_mono, which(STAR_cells_in_gene<3454, arr.ind=TRUE))
  sorted_distance_mono_highlycovered <- sorted_distance_with_mono[-highlycovered_mono, ]
  
  lesscovered_nomono <- with(sorted_distance_no_mono, which((STAR_cells_in_gene<=3454 && STAR_cells_in_gene>=1381) , arr.ind=TRUE))
  print(lesscovered_nomono)
  sorted_distance_nomono_lesscovered <- sorted_distance_no_mono[lesscovered_nomono, ]
  
  print(sorted_distance_nomono_highlycovered[1:15,]) 
  print(sorted_distance_mono_highlycovered[1:15,])
  print(sorted_distance_nomono_lesscovered[1:100,])
  
  # s1<- sorted1[1:100,]
  # #s2<- sorted2[1:500,]
  # s3<- sorted3[1:500,]
  # 
  # final_candidates <- Reduce(intersect, list(s2,s3))
  
  
  # final_candidates = list()
  # print(final_candidates)
  # x <- 1
  # print(length(final_candidates))
  # while (x < 51) {
  #   
  #   for (i in 1:nrow(s1)){
  # 
  #     for (j in 1:nrow(s2)){
  # 
  #       if (s1$ms_start[i] == s2$ms_start[j] && s1$ms_end[i]== s2$ms_end[j] && s1$specific_gene_ncbi_id[i] == s2$specific_gene_ncbi_id[j])
  #         {
  #           final_candidates[x] <-s2$ms_start[j]
  #         }
  #     }
  #     
  #     x = x+1 
  #   }
  # }
  # 
  #print(final_candidates)
}

find_candidates("/Users/krystenharvey/Documents/Landau Lab/Lineage_Project/new_candidate_MS_loci - new_candidate_MS_loci.csv")