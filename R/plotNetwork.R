plotNetwork <- function(fit, threshold=0.05, thresh.type="p.val", layout="circular", legend.pos="left") {
#  Input:
####  - fit : a complete DINGO or iDINGO result
####  - threshold : threshold (either a p-value or absolute differential score) for inclusion in differential network
####  - thresh.type : either "p.val" for p-value, or "diff.score" for absolute differential score
####  - layout : Either "circular" (which automatically becomes cylindrical when multi-platform) or one of
####    igraph's supported layouts
  ##  - legend.pos: Position of legend, in c("left", "right", "top", "bottom")

if (!(thresh.type %in% c("p.val", "diff.score"))) {
  stop("thresh.type must be either 'p.val' or 'diff.score'")
}

fit.type <- ifelse(length(fit)==6, yes = "iDINGO", no = "DINGO")

data.for.plot <- cbind(fit$genepair, fit$diff.score, fit$p.val,
                         ifelse(fit$R1 * fit$R2 > 0, yes=1, no=0))
colnames(data.for.plot)[3:5] <- c("score", "p.val", "conserved")

if (thresh.type == "p.val") {
  data.for.plot <- data.for.plot[data.for.plot$p.val < threshold,]
} else {
  data.for.plot <- data.for.plot[abs(data.for.plot$score) > threshold,]
}

rownames(data.for.plot) <- NULL
network.for.plot <- graph.data.frame(data.for.plot, directed=FALSE)

E(network.for.plot)$color <- ifelse(E(network.for.plot)$score>0, "red", "blue") #red if positive, blue if negative
E(network.for.plot)$width <- 0.4*(abs(E(network.for.plot)$score)-1.9)
E(network.for.plot)$lty <- E(network.for.plot)$conserved

V(network.for.plot)$size <- 40*degree(network.for.plot) / max(degree(network.for.plot))
V(network.for.plot)$size[V(network.for.plot)$size < 20] <- 20

# This bit from http://kieranhealy.org/blog/archives/2011/02/18/aligning-labels-in-circular-igraph-layouts/
# For circular layout
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(rescale(x, c(0, 2 * pi), range(x)))
}
cylinder.layout <- function(radii, angles, z) {
  cbind(radii*cos(angles)+z, -radii*sin(angles)-z/16)
}

# For v1.0.2 - $name resassignment can cause an error if 2 elements from
# different platforms have the same name. Commenting out...
if (fit.type == "iDINGO") {
  V(network.for.plot)$group <- gsub("\\|\\|.*", "", V(network.for.plot)$name)
#  V(network.for.plot)$name <- gsub(".*\\|\\|", "", V(network.for.plot)$name)
  groups <- unique(V(network.for.plot)$group)
}


network.for.plot <- permute(network.for.plot, rank(V(network.for.plot)$name))

if (fit.type == "iDINGO") {
  angles <- rep.int(0, vcount(network.for.plot))
  plat_a <- which(V(network.for.plot)$group == groups[1])
  plat_b <- which(V(network.for.plot)$group == groups[2])
  angles[plat_a] <- (1:length(plat_a)-1) * 2 * pi / length(plat_a)
  angles[plat_b] <- (1:length(plat_b)-1) * 2 * pi / length(plat_b)
  z <- ifelse(V(network.for.plot)$group == groups[1], -3, 3)
  z[V(network.for.plot)$group == groups[2]] <- 0

  if (length(groups) == 3) {
    plat_c <- which(V(network.for.plot)$group == groups[3])
    angles[plat_c] <- (1:length(plat_c)-1) * 2 * pi / length(plat_c)
  }

  radii <- 1
  layout.n <- cylinder.layout(radii, angles, z)
} else { layout.n <- layout_in_circle(network.for.plot) }




if (layout == "circular") {
  visNet <- visIgraph(network.for.plot, layout = "layout.norm", layoutMatrix = layout.n, physics = FALSE)
} else {
  visNet <- visIgraph(network.for.plot, layout = layout, physics = FALSE)
}

if (fit.type == "iDINGO") {
  visNet <- visNet %>%
    visGroups(groupname = groups[1], color = list(background = "lightgreen", border = "black", highlight = "black")) %>%
    visGroups(groupname = groups[2], color = list(background = "orange", border = "black", highlight = "black"))
  if (length(groups) == 3) {
    visNet <- visNet %>%
      visGroups(groupname = groups[3], color = list(background = "violet", border = "black", highlight = "black"))
  }
}

visNet %>%
  visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE),
                        nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = FALSE) %>%
  visLegend(position = legend.pos)

}
