options(stringsAsFactors =  F)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(stringdist)
library(pals)
# set working directory to the directory that contains the input files.
# setwd("Path/to/folder/")


#experiment number, script is looking for e. g. "./exp_1.csv" if exp is 1
exp = 1
# toggle wheter or not the grouping should be done with graphs.
# with graphs its easier and delivers better results.
use.graphs = F

# defining label colours
my.colors = c(cols25(), alphabet(), alphabet2())

# reading data.
infilename = paste0("./exp_", exp, ".csv")
print(paste("reading from", infilename))
read.data = read.csv(infilename)

#column number that holds the first quant value.
quant.start.index = which(colnames(read.data) == "q126")
#column number that holds the last quant value
quant.end.index = which(colnames(read.data) == "q131")

# remove rows that have only NA quant values.
read.data = read.data[apply(read.data[,quant.start.index:quant.end.index], 1, function(x){!(all(is.na(x)))}), ]
# make all sequences uppercase and replace I and J with L
read.data$sequence = gsub("(L|I)", "J", read.data$sequence)

################################################################################
# optional: read sample names, which can then be displayed for each heatmap.
read.sample.names = read.delim("./cellline_names.txt")
sample.names = read.sample.names[exp, 2:11]
################################################################################

# creating list of proteins included in the current file.
protein.ids = sort(unique(unlist(lapply(read.data$proteinID, function(x) {strsplit(x, "/")}))))

################################################################################
# optional: read protein ids from different file and intersect them with the current list
# this reduces the number of heatmaps in the output and allows for filtering 
read.protein.ids = read.delim("./reduced_protein_ids.txt")
protein.ids = intersect(protein.ids, read.protein.ids$Entry)
################################################################################

################################################################################
# optional: read possibly existing data containing information on already identified variants. 
# this will put VAR_ in front of the PSMs that are in the input file and has a hamming distance of 1 or less
# compared to a variant found in the variant file. the variant file consists of just the variant peptides. 
variants = read.delim("./var_peps.txt", header = F)
variants$upperSeq = gsub("(L|I)", "J", variants$V1)

has.potential.variant.seq = unlist(lapply(read.data$sequence, function(x) {
  any(stringdist(x, variants$upperSeq, method = "hamm") <= 1)
}))

read.data$pot.var.seq = has.potential.variant.seq
################################################################################

# start & end index in the sorted protein ids 
kmin = 1
kmax = length(protein.ids)

# normalisation
quant.values = read.data[,quant.start.index:quant.end.index]
norm.quant.values = t(apply(quant.values, 1, function(x){
  x/max(x, na.rm = T)
}))
read.data[,quant.start.index:quant.end.index] = norm.quant.values

# file writer
outfilename = paste0("./out_", ifelse(use.graphs, "graphs", "conv"), "_exp_", exp, ".pdf")
print(paste("writing to", outfilename))
pdf(outfilename,  width=30, height=40, onefile = T)

# main loop, for each protein
for(k in kmin:kmax) {
  #progress report
  if ((k%%10==0 )) {print(paste0(k, "/", kmax, " (", round(k/kmax, 4)*100, "%)"))}
  # subsetting the data, taking into account only psms that map to the specified protein.
  grepd.data = read.data[grep(protein.ids[k], read.data$proteinID),]
  
  #filtering: to few PSMs make analysis impossible, too many make the heatmaps unreadable.
  if (nrow(grepd.data) > 2 && nrow(grepd.data) < 50){
    # ordering data by psm id
    grepd.data = grepd.data[order(grepd.data$psmID), ]
    
    #split SMiV and identification PSMs
    starts.with.PSM = startsWith(grepd.data$psmID, "PSM")
    ident.length = length(which(!starts.with.PSM))
    smiv.length = length(which(starts.with.PSM))
    
    #filter: if there is anything to compare in the data.
    if (ident.length > 1 && smiv.length > 0 ) {
      # create new (shorter) names to display in the heatmaps
      new.names.ident = data.frame(psmID = grepd.data[!starts.with.PSM, "psmID"], newPsmID = paste0("P", 1:ident.length )) 
      new.names.smiv = data.frame(psmID = grepd.data[starts.with.PSM, "psmID"], newPsmID = paste0("S", 1:smiv.length ))
      new.names = rbind(new.names.ident, new.names.smiv)
      grepd.data$newPsmID = paste0(ifelse(grepd.data$pot.var.seq, "VAR_", ""), new.names$newPsmID)
      
      if (use.graphs) {
        library(igraph)
        # find subgroups that belong together according to SMiV
        # find which original psms the smiv ids map to 
        occurence.in.mapping = sapply(grepd.data[, "refpsmID"], function(x) {
          match(new.names.ident$psmID, unlist(strsplit(x, "/")))
        }, USE.NAMES = F)
        
        # find all peptides that share a sequence (stringdistance <=1) (optimal string alignment)
        stringdistances = sapply(grepd.data[, "sequence"], function(x){
          stringdist(x, grepd.data[!starts.with.PSM, "sequence"])
        }, USE.NAMES = F)
        
        #combining the information extracted above.
        combined.matrix = (stringdistances<=1) | !is.na(occurence.in.mapping)
        colnames(combined.matrix) = new.names$newPsmID
        rownames(combined.matrix) = new.names.ident$newPsmID
        
        # transform data to long format
        melted.matrix = melt(combined.matrix)
        #remove the last column
        combinations = melted.matrix[melted.matrix$value, -3]
        
        #create network and extract groups from there.
        network = graph_from_data_frame(combinations, F)
        membership = components(network)$membership
        groups = data.frame(newPsmID = names(membership), group = as.integer(membership))
        
      } else {
        #original method of calculating membership, relatively complicated. 
        smiv.identifications = grepd.data[starts.with.PSM, c(2, 7)]
        other.identifications = grepd.data[!starts.with.PSM, 2]
        leftover.identifications = other.identifications
        
        
        combinations = data.frame(i = integer(),smiv = character(), other = character())
        index = 1
        for (x in 1:nrow(smiv.identifications)){
          split = unlist(strsplit(smiv.identifications[x, "refpsmID"], "/"))
          in.both = intersect(split, other.identifications)
          if (length(in.both) == 0) {
            combinations[index, ] = c(index, smiv.identifications[x, "psmID"], paste0("NOMATCH", index))
            index = index+1
          }else{
            for (y in in.both) {
              combinations[index, ] = c(index, smiv.identifications[x, "psmID"], y)
              index = index + 1
              leftover.identifications = setdiff(leftover.identifications, split)
            }
          }
        }
        for (x in leftover.identifications) {
          combinations[index, ] = c(index, NA, x)
          index = index + 1
        }
        
        grouped = data.frame(group=integer, psmID = character())
        group.id = 1
        smiv.ids = sort(unique(combinations[,2]))
        for (smiv in smiv.ids) {
          curr = combinations[is.element(combinations$smiv, smiv),]
          initial.rows = nrow(curr)
          while(nrow(curr) > 0){
            to.remove = curr$i
            curr.psms = unlist(c(as.matrix(curr[,c("smiv", "other")])))
            grouped = rbind(grouped, data.frame(group = rep(group.id, length(curr.psms)),  psmID = curr.psms))
            
            other.ids = unique(curr$other)
            combinations = combinations[!is.element( combinations$i, to.remove), ]
            curr = combinations[is.element(combinations$other, other.ids),]
            curr.psms = unlist(c(as.matrix(curr[,2:3])))
            grouped = rbind(grouped, data.frame(group = rep(group.id, length(curr.psms)),  psmID = curr.psms))
            
            to.remove = curr$i
            smiv = curr$smiv
            combinations = combinations[!is.element(combinations$i, to.remove), ]
            curr = combinations[is.element(combinations$smiv, smiv),]
          }
          if (initial.rows > 0) {
            group.id = group.id + 1
          }
        }
        
        if (nrow(combinations) > 0) {
          for (others in 1:nrow(combinations)){
            grouped = rbind(grouped, data.frame(group = group.id, psmID = combinations[others, "other"]))
            group.id = group.id + 1
          }
        }
        grouped = unique(grouped[!startsWith(grouped[,2], "NOMATCH"),])
        groups = data.frame(newPsmID = new.names[match(grouped$psmID, new.names$psmID), "newPsmID"], group = grouped$group)
      }
      
      #add the group to the data and reorder columns. 
      grepd.data = merge(grepd.data, groups, all.x = T)
      if (any(is.na(grepd.data$group))){
        ungrouped.elements = length(grepd.data[is.na(grepd.data$group), "group"])
        first.new.group = (max(grepd.data$group, na.rm = T)+1)
        grepd.data[is.na(grepd.data$group), "group"] = first.new.group:(first.new.group + ungrouped.elements -1)
      }
      grepd.data = grepd.data[order(grepd.data$psmID), c(2:20, 1, 21)]
      grepd.data = grepd.data[order(grepd.data$group, grepd.data$psmID), ]
      
      # extract just quant values.
      grepd.data.matrix = grepd.data[,quant.start.index:quant.end.index]
      
      nelms = nrow(grepd.data)
      plots = vector("list", 10)
      colors = my.colors[grepd.data$group]
      
      #inner loop for each cell line.
      for(col.id in 1:ncol(grepd.data.matrix)) {
        # extract one column
        current.sample.data = grepd.data.matrix[,col.id]
        #create the matrix to fill
        ratios = matrix(numeric(nelms*nelms), ncol = nelms)
        rownames(ratios) = colnames(ratios) = grepd.data$newPsmID
        
        #ratio calculation, heart of the approach.
        for (x in 1:nelms) {
          for (y in 1:nelms){
            ratios[x, y] = log2(current.sample.data[x]/current.sample.data[y])
          }
        }
        
        #transform data for plotting
        melted.ratios = melt(ratios)
        melted.ratios[,1] = as.factor(melted.ratios[,1])
        melted.ratios[,2] = as.factor(melted.ratios[,2])
        
        #save all 10 sample plots.
        plots[[col.id]] = 
          ggplot(data = melted.ratios, aes(x = Var1, y = Var2, fill = value)) + 
          geom_tile()  +
          scale_fill_gradient2(low="red", high = "blue", mid ="white", midpoint = 0, na.value = "darkgrey")  + 
          theme(legend.position="none") + 
          theme(axis.title = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
          theme(axis.text = element_text(color=colors, face = "bold", size = 16)) + 
          ggtitle(paste(sample.names[col.id])) + 
          theme(title = element_text(size = 18, face = "bold"))
      }
      
      #print all plots to one page.
      grid.arrange(grobs = plots, ncol = 3, top=textGrob(protein.ids[k], gp=gpar(fontsize=30, fontface="bold")))
    }
  }
}
dev.off()
