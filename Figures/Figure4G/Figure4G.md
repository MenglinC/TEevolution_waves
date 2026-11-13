# G

## G.1
```
#MA1630.3

```


## G.2.1

```shell
streme --dna --p /mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/results/motif/SRR3112448_summits.top1000.flank100.fa --oc /mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/results/motif_0806/ --minw 6 --maxw 10 --nmotifs 10
```

## G.2.2

```r

**############################################################
# Plot STREME motif logos
# Forward + Reverse complement
# One-page PDF
############################################################


library(seqLogo)
library(grid)
library(gridExtra)


############################################################
# Input
############################################################

streme_file <- "streme.txt"

output_pdf <- "STREME_all_motifs_logo.pdf"



############################################################
# Read STREME motif matrices
############################################################

lines <- readLines(streme_file)


motif_start <- grep("^MOTIF", lines)


motifs <- list()
motif_info <- list()


for(i in seq_along(motif_start)){


  start <- motif_start[i]


  # matrix line
  mat_line <- grep(
    "letter-probability matrix",
    lines[(start+1):length(lines)]
  )[1] + start



  # width
  w <- as.numeric(
    sub(
      ".*w= ([0-9]+).*",
      "\\1",
      lines[mat_line]
    )
  )


  # sites
  sites <- as.numeric(
    sub(
      ".*nsites= ([0-9]+).*",
      "\\1",
      lines[mat_line]
    )
  )


  # E value
  evalue <- sub(
    ".*E= ([^ ]+).*",
    "\\1",
    lines[mat_line]
  )



  # read matrix

  mat <- do.call(
    rbind,
    lapply(
      lines[(mat_line+1):(mat_line+w)],
      function(x){
        as.numeric(
          strsplit(trimws(x),"\\s+")[[1]]
        )
      }
    )
  )


  colnames(mat) <- c("A","C","G","T")


  motifs[[i]] <- mat


  motif_info[[i]] <- list(
    width=w,
    sites=sites,
    E=evalue
  )

}



############################################################
# Reverse complement
############################################################


revcomp_matrix <- function(mat){

  rc <- mat[,c("T","G","C","A")]

  rc <- rc[nrow(rc):1,]

  colnames(rc) <- c("A","C","G","T")

  return(rc)

}




############################################################
# Generate logo grobs
############################################################


logo_list <- list()


for(i in seq_along(motifs)){


  mat <- motifs[[i]]


  #################################
  # Forward
  #################################

  pwm_forward <- makePWM(t(mat))


  p1 <- grid.grabExpr({

    seqLogo(
      pwm_forward,
      ic.scale=TRUE
    )

  })


  title1 <- textGrob(
    paste0(
      "Motif ",
      i,
      " Forward\n",
      "w=",
      motif_info[[i]]$width,
      "  sites=",
      motif_info[[i]]$sites,
      "  E=",
      motif_info[[i]]$E
    ),
    gp=gpar(fontsize=10)
  )



  #################################
  # Reverse complement
  #################################


  mat_rc <- revcomp_matrix(mat)


  pwm_rc <- makePWM(t(mat_rc))


  p2 <- grid.grabExpr({

    seqLogo(
      pwm_rc,
      ic.scale=TRUE
    )

  })


  title2 <- textGrob(
    paste0(
      "Motif ",
      i,
      " Reverse complement"
    ),
    gp=gpar(fontsize=10)
  )



  logo_list[[length(logo_list)+1]] <-
    arrangeGrob(
      title1,
      p1,
      heights=c(0.15,1)
    )


  logo_list[[length(logo_list)+1]] <-
    arrangeGrob(
      title2,
      p2,
      heights=c(0.15,1)
    )

}



############################################################
# Plot one page
############################################################


pdf(
  output_pdf,
  width=12,
  height=10
)


grid.arrange(
  grobs=logo_list,
  ncol=2
)


dev.off()



cat(
  "Finished:",
  output_pdf,
  "\n"
)**
```
