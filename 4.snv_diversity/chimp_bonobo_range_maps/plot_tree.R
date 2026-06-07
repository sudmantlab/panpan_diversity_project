setwd("/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN")
# Load required packages
library(ggplot2)
library(ape)
library(ggtree)
#Western-NC: 250,000 years (marked as δ)
#Central-Eastern: 139,000-182,000 years (median = 160,500 years)
#All chimps coalesce: 544,000-633,000 years (median = 588,500 years)
#Chimp-Bonobo split: 1,660,000-2,100,000 years (median = 1,880,000 years)
#Human split: 6,000,000 years

#Newick format:
#(Human:6000000,(Bonobo:1880000,((Western:250000,NC:250000):338500,(Central:160500,Eastern:160500):428000):1291500):4120000);


newick_str <- "(Human:6000000,(Bonobo:1880000,((Western:250000,NC:250000):338500,(Central:160500,Eastern:160500):428000):1291500):4120000);"

# Read the tree
tree <- read.tree(text = newick_str)

# Create curved tree plot
p <- ggtree(tree, layout="rectangular", branch.length='length') +
  # Add curved branches
  geom_tree(layout="rectangular", 
            lineend="round",    # round line endings
            linejoin="round",   # round line joins
            size=3) +           # thicker lines
  # Add tip labels
  geom_tiplab(aes(color=label), 
              size=5,           # larger text
              fontface="bold",  # bold text
              hjust=0.5) +
  # Set custom colors
  scale_color_manual(values=c(
    "Western" = "#9dced9",
    "NC" = "#fe604c",
    "Eastern" = "#ffb35a",
    "Central" = "#4c5d4c",
    "Bonobo" = "salmon",       
    "Human" = "black")) +
  # Rotate and flip
  coord_flip() + 
  scale_x_reverse() +
  # Remove legend
  theme(legend.position="none") +
  # Add time scale
  scale_y_continuous(
    name="Time (kya)",
    breaks=seq(0, 2000000, by=500000),
    labels=function(x) paste0(x/1000, " kya"),
    sec.axis = dup_axis(name = "Time (Mya)",
                        labels = function(x) paste0(x/1000000, " Mya"))
  )

# Print the plot
print(p) 