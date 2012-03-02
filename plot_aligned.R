source("data_handling.R")

# Plots aligned sequences numbers seqnum1 and seqnum2 into a pdf
# file. (The required aligned data files must be in the aligneddata
# directory.)
plot_aligned <- function(seqnum1, seqnum2, outputfile = "test.pdf") {
    pdf(file = outputfile, height = 5, width = 7)
    par(mfcol = c(2, 1), mar=c(5, 2, 2, 1))

    draw_plot(seqnum1, seqnum2, sprintf("Sequence %d", seqnum1))
    draw_plot(seqnum2, seqnum1, sprintf("Sequence %d", seqnum2))

    box()
    dev.off()
}

draw_plot <- function(seqnum1, seqnum2, mainlab) {
    sequence <- mean_normalize(load_aligned(seqnum1, seqnum2))
    colors = rainbow(ncol(sequence))

    # Bug(?): for some reason ylab does not appear in the resulting pdf.
    plot(sequence[,1], type = "l", col = colors[1], ylim = c(min(sequence), max(sequence)),
        xlim = c(0, nrow(sequence)), xlab = "Frames", ylab = "Features",
        main = mainlab)
    for (i in 2:ncol(sequence))
        lines(sequence[,i], col = colors[i])
}

mean_normalize <- function(sequence) {
    return(t(apply(sequence, 1, '-', apply(sequence, 2, mean))))
}
