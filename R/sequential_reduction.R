#  ----------------------------------------------------------------------
# Figure 1: sequential reduction as vertex elimination on the variable
# graph of a five-variable chain.
#
#   (a) the target sum and the graph it induces
#   (b) interior-first elimination, which creates fill-in
#   (c) endpoint-first elimination, which does not
#   (d) operation counts against brute-force enumeration
#
# Base R only. Drawing happens in a fixed pixel coordinate system with the
# y axis reversed, so every coordinate below reads top-down like the page.

# ---- geometry -----------------------------------------------------------

fig_width_px <- 1120
fig_height_px <- 950

panel_margin <- 40
panel_gap <- 34
panel_width <- (fig_width_px - 2 * panel_margin - 2 * panel_gap) / 3
panel_left <- panel_margin + (0:2) * (panel_width + panel_gap)
node_spacing <- 52
node_offset <- 69
node_radius <- 13

y_setup_node <- 128
y_setup_eq <- 88
y_setup_note <- 176

y_b_eq <- 284
y_b_node <- 348
y_b_caption <- 430

y_c_eq <- 538
y_c_node <- 602
y_c_caption <- 684

plot_left <- 560
plot_right <- 1080
plot_top <- 742
plot_bottom <- 902

log_min <- 2
log_max <- 10
n_min <- 3
n_max <- 10
k_states <- 10

# ---- palette ------------------------------------------------------------

col_ink <- "#111111"
col_fillin <- "#3f3f3f"
col_halo <- "#cbcbcb"
col_gone_line <- "#bcbcbc"
col_gone_text <- "#9a9a9a"
col_caption <- "#3a3a3a"
col_rule <- "#dcdcdc"
col_purple <- "#8E44AD" # RTMB accent, matches models/open_nmixture_spde.R

# ---- unit helpers -------------------------------------------------------

# One user unit is one point on the svg device, and lwd 1 is 1/96 inch.
lwd_from_pt <- function(width_pt) {
  value <- width_pt / 0.75
  value
}

# Base pointsize is set to 12, so cex is the requested size over 12.
cex_from_pt <- function(size_pt) {
  value <- size_pt / 12
  value
}

node_x <- function(panel, index) {
  base <- panel_left[panel]
  value <- base + node_offset + (index - 1) * node_spacing
  value
}

panel_centre <- function(panel) {
  value <- panel_left[panel] + panel_width / 2
  value
}

# ---- expression builders ------------------------------------------------

pair_factor <- function(letter, subscript, first, second) {
  head_sym <- as.name(letter)
  expr <- bquote(paste(
    .(head_sym)[.(subscript)], "(", x[.(first)], ", ", x[.(second)], ")"
  ))
  expr
}

single_factor <- function(letter, subscript, first) {
  head_sym <- as.name(letter)
  expr <- bquote(paste(.(head_sym)[.(subscript)], "(", x[.(first)], ")"))
  expr
}

join_two <- function(left_expr, right_expr) {
  expr <- bquote(paste(.(left_expr), "  ", .(right_expr)))
  expr
}

with_equals <- function(left_expr) {
  expr <- bquote(paste(.(left_expr), "  =  "))
  expr
}

# ---- primitive drawing --------------------------------------------------

circle_theta <- seq(0, 2 * pi, length.out = 96)

draw_disc <- function(x_centre, y_centre, fill, border, dash, width_pt) {
  x_ring <- x_centre + node_radius * cos(circle_theta)
  y_ring <- y_centre + node_radius * sin(circle_theta)
  line_width <- lwd_from_pt(width_pt)
  polygon(x_ring, y_ring,
    col = fill, border = border, lwd = line_width,
    lty = dash
  )
  invisible(NULL)
}

draw_straight <- function(x_start, x_end, y_line, colour, width_pt, dash) {
  line_width <- lwd_from_pt(width_pt)
  segments(x_start, y_line, x_end, y_line,
    col = colour, lwd = line_width,
    lty = dash, lend = "round"
  )
  invisible(NULL)
}

# Quadratic Bezier with the control point pushed below the node line, so
# fill-in edges nest instead of overprinting the chain.
arc_points <- function(x_start, x_end, y_line, span) {
  depth <- 26 * span
  x_control <- (x_start + x_end) / 2
  y_control <- y_line + depth
  t_grid <- seq(0, 1, length.out = 80)
  weight_a <- (1 - t_grid)^2
  weight_b <- 2 * (1 - t_grid) * t_grid
  weight_c <- t_grid^2
  x_curve <- weight_a * x_start + weight_b * x_control + weight_c * x_end
  y_curve <- weight_a * y_line + weight_b * y_control + weight_c * y_line
  result <- list(x = x_curve, y = y_curve)
  result
}

draw_arc <- function(x_start, x_end, y_line, span, colour, width_pt, dash) {
  curve <- arc_points(x_start, x_end, y_line, span)
  line_width <- lwd_from_pt(width_pt)
  lines(curve$x, curve$y,
    col = colour, lwd = line_width, lty = dash,
    lend = "round"
  )
  invisible(NULL)
}

draw_chevron <- function(x_left, y_centre) {
  x_path <- c(x_left, x_left + 7, x_left)
  y_path <- c(y_centre - 6, y_centre, y_centre + 6)
  line_width <- lwd_from_pt(1.6)
  lines(x_path, y_path, col = "#a8a8a8", lwd = line_width, lend = "round")
  invisible(NULL)
}

draw_rule <- function(y_line) {
  line_width <- lwd_from_pt(1)
  segments(40, y_line, 1080, y_line, col = col_rule, lwd = line_width)
  invisible(NULL)
}

# ---- text ---------------------------------------------------------------

# y is the text baseline, matching the layout constants above.
draw_text <- function(x_at, y_at, label, size_pt, colour, align,
                      font_family, rotation = 0) {
  size_cex <- cex_from_pt(size_pt)
  text(x_at, y_at,
    labels = label, cex = size_cex, col = colour,
    adj = c(align, 0), family = font_family, srt = rotation
  )
  invisible(NULL)
}

draw_prose <- function(x_at, y_at, label, size_pt, colour, align = 0,
                       rotation = 0) {
  draw_text(x_at, y_at, label, size_pt, colour, align, "sans", rotation)
  invisible(NULL)
}

draw_math <- function(x_at, y_at, label, size_pt, colour, align = 0.5,
                      font_family = "serif") {
  draw_text(x_at, y_at, label, size_pt, colour, align, font_family)
  invisible(NULL)
}

math_width <- function(label, size_pt, font_family = "serif") {
  size_cex <- cex_from_pt(size_pt)
  value <- strwidth(label, cex = size_cex, family = font_family)
  value
}

# An equation is laid out as four runs so the summation sign can carry an
# index set below and to the right of it, the way it appears in the text.
draw_equation <- function(x_centre, y_baseline, lhs_expr, index_expr,
                          rhs_expr, size_pt) {
  sigma_expr <- quote(Sigma)
  size_sigma <- size_pt * 1.35
  size_index <- size_pt * 0.68
  width_lhs <- math_width(lhs_expr, size_pt)
  width_sigma <- math_width(sigma_expr, size_sigma)
  width_index <- math_width(index_expr, size_index)
  width_rhs <- math_width(rhs_expr, size_pt)
  gap <- 3
  total <- width_lhs + width_sigma + width_index + gap + width_rhs
  x_cursor <- x_centre - total / 2
  draw_math(x_cursor, y_baseline, lhs_expr, size_pt, col_ink, 0)
  x_cursor <- x_cursor + width_lhs
  draw_math(x_cursor, y_baseline + 1, sigma_expr, size_sigma, col_ink, 0)
  x_cursor <- x_cursor + width_sigma
  y_index <- y_baseline + size_pt * 0.74
  draw_math(x_cursor, y_index, index_expr, size_index, col_ink, 0)
  x_cursor <- x_cursor + width_index + gap
  draw_math(x_cursor, y_baseline, rhs_expr, size_pt, col_ink, 0)
  invisible(NULL)
}

# ---- graph frames -------------------------------------------------------

chain_edges <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 5))

draw_graph <- function(panel, frame, y_line, accent = col_ink) {
  for (item in frame$halo) {
    x_start <- node_x(panel, item$from)
    x_end <- node_x(panel, item$to)
    span <- item$to - item$from
    if (identical(item$kind, "straight")) {
      draw_straight(x_start, x_end, y_line, col_halo, 10, "solid")
    } else {
      draw_arc(x_start, x_end, y_line, span, col_halo, 10, "solid")
    }
  }
  for (item in frame$solid) {
    x_start <- node_x(panel, item[1])
    x_end <- node_x(panel, item[2])
    draw_straight(x_start, x_end, y_line, col_ink, 1.7, "solid")
  }
  for (item in frame$fill) {
    x_start <- node_x(panel, item[1])
    x_end <- node_x(panel, item[2])
    span <- item[2] - item[1]
    draw_arc(x_start, x_end, y_line, span, col_fillin, 1.5, "32")
  }
  for (index in 1:5) {
    x_node <- node_x(panel, index)
    is_active <- isTRUE(frame$active == index)
    is_gone <- index %in% frame$gone
    if (is_active) {
      draw_disc(x_node, y_line, accent, col_ink, "solid", 1.6)
      label_colour <- "#ffffff"
    } else if (is_gone) {
      draw_disc(x_node, y_line, "#ffffff", col_gone_line, "33", 1.3)
      label_colour <- col_gone_text
    } else {
      draw_disc(x_node, y_line, "#ffffff", col_ink, "solid", 1.6)
      label_colour <- col_ink
    }
    node_label <- bquote(x[.(index)])
    draw_math(x_node, y_line + 4.5, node_label, 14, label_colour, 0.5)
  }
  invisible(NULL)
}

halo_item <- function(kind, from, to) {
  value <- list(kind = kind, from = from, to = to)
  value
}

# ---- panel content ------------------------------------------------------

frame_setup <- list(
  active = NA,
  gone = integer(0),
  solid = chain_edges,
  fill = list(),
  halo = list()
)

b_lhs_1 <- with_equals(pair_factor("g", 13, 1, 3))
b_rhs_1 <- join_two(pair_factor("f", 12, 1, 2), pair_factor("f", 23, 2, 3))
b_lhs_2 <- with_equals(pair_factor("g", 35, 3, 5))
b_rhs_2 <- join_two(pair_factor("f", 34, 3, 4), pair_factor("f", 45, 4, 5))
b_lhs_3 <- with_equals(pair_factor("g", 15, 1, 5))
b_rhs_3 <- join_two(pair_factor("g", 13, 1, 3), pair_factor("g", 35, 3, 5))

c_lhs_1 <- with_equals(single_factor("h", 2, 2))
c_rhs_1 <- pair_factor("f", 12, 1, 2)
c_lhs_2 <- with_equals(single_factor("h", 3, 3))
c_rhs_2 <- join_two(single_factor("h", 2, 2), pair_factor("f", 23, 2, 3))
c_lhs_3 <- with_equals(single_factor("h", 4, 4))
c_rhs_3 <- join_two(single_factor("h", 3, 3), pair_factor("f", 34, 3, 4))

note_b_1 <- bquote(paste(
  K^3, " operations, clique {", x[1], ", ", x[2],
  ", ", x[3], "}"
))
note_b_2 <- bquote(paste(
  K^3, " operations, clique {", x[3], ", ", x[4],
  ", ", x[5], "}"
))
note_b_3 <- bquote(paste(
  K^3, " operations, clique {", x[1], ", ", x[3],
  ", ", x[5], "}"
))
note_c_1 <- bquote(paste(
  K^2, " operations, clique {", x[1], ", ", x[2],
  "}"
))
note_c_2 <- bquote(paste(
  K^2, " operations, clique {", x[2], ", ", x[3],
  "}"
))
note_c_3 <- bquote(paste(
  K^2, " operations, clique {", x[3], ", ", x[4],
  "}"
))

cap_b_1 <- bquote(paste("fill-in edge ", x[1], "-", x[3]))
cap_b_2 <- bquote(paste("fill-in edge ", x[3], "-", x[5]))
cap_b_3 <- bquote(paste("fill-in edge ", x[1], "-", x[5]))
cap_c_1 <- bquote(paste(
  "no fill-in: ", h[2], " depends on ", x[2],
  " alone"
))
cap_c_2 <- bquote(paste(
  "no fill-in: ", h[3], " depends on ", x[3],
  " alone"
))
cap_c_3 <- bquote(paste(
  "no fill-in: ", h[4], " depends on ", x[4],
  " alone"
))

row_b <- list(
  list(
    active = 2,
    gone = integer(0),
    solid = chain_edges,
    fill = list(c(1, 3)),
    halo = list(
      halo_item("straight", 1, 2), halo_item("straight", 2, 3),
      halo_item("arc", 1, 3)
    ),
    lhs = b_lhs_1, index = quote(x[2]), rhs = b_rhs_1,
    note = note_b_1, caption = cap_b_1
  ),
  list(
    active = 4,
    gone = c(2),
    solid = list(c(3, 4), c(4, 5)),
    fill = list(c(1, 3), c(3, 5)),
    halo = list(
      halo_item("straight", 3, 4), halo_item("straight", 4, 5),
      halo_item("arc", 3, 5)
    ),
    lhs = b_lhs_2, index = quote(x[4]), rhs = b_rhs_2,
    note = note_b_2, caption = cap_b_2
  ),
  list(
    active = 3,
    gone = c(2, 4),
    solid = list(),
    fill = list(c(1, 3), c(3, 5), c(1, 5)),
    halo = list(
      halo_item("arc", 1, 3), halo_item("arc", 3, 5),
      halo_item("arc", 1, 5)
    ),
    lhs = b_lhs_3, index = quote(x[3]), rhs = b_rhs_3,
    note = note_b_3, caption = cap_b_3
  )
)

row_c <- list(
  list(
    active = 1,
    gone = integer(0),
    solid = chain_edges,
    fill = list(),
    halo = list(halo_item("straight", 1, 2)),
    lhs = c_lhs_1, index = quote(x[1]), rhs = c_rhs_1,
    note = note_c_1, caption = cap_c_1
  ),
  list(
    active = 2,
    gone = c(1),
    solid = list(c(2, 3), c(3, 4), c(4, 5)),
    fill = list(),
    halo = list(halo_item("straight", 2, 3)),
    lhs = c_lhs_2, index = quote(x[2]), rhs = c_rhs_2,
    note = note_c_2, caption = cap_c_2
  ),
  list(
    active = 3,
    gone = c(1, 2),
    solid = list(c(3, 4), c(4, 5)),
    fill = list(),
    halo = list(halo_item("straight", 3, 4)),
    lhs = c_lhs_3, index = quote(x[3]), rhs = c_rhs_3,
    note = note_c_3, caption = cap_c_3
  )
)

draw_row <- function(row, y_eq, y_node, y_caption, accent = col_ink) {
  for (panel in seq_along(row)) {
    frame <- row[[panel]]
    x_centre <- panel_centre(panel)
    draw_equation(x_centre, y_eq, frame$lhs, frame$index, frame$rhs, 14.5)
    draw_math(
      x_centre, y_eq + 28, frame$note, 12, col_gone_text,
      0.5, "sans"
    )
    draw_graph(panel, frame, y_node, accent)
    draw_math(
      x_centre, y_caption, frame$caption, 12.5, col_caption,
      0.5, "sans"
    )
    if (panel < 3) {
      draw_chevron(panel_left[panel] + panel_width + panel_gap / 2 - 3, y_node)
    }
  }
  invisible(NULL)
}

# ---- operation counts ---------------------------------------------------

ops_brute <- function(n, k) {
  value <- (n - 1) * k^n
  value
}

ops_interior <- function(n, k) {
  value <- (n - 2) * k^3
  value
}

ops_endpoint <- function(n, k) {
  value <- (n - 1) * k^2
  value
}

count_rows <- list(
  list(
    name = "brute force", dash = "solid", width_pt = 2,
    formula = quote(phantom() %~~% paste("(", n - 1, ")", K^n)),
    fun = ops_brute
  ),
  list(
    name = "interior-first", dash = "43", width_pt = 1.6,
    formula = quote(phantom() %~~% paste("(", n - 2, ")", K^3)),
    fun = ops_interior
  ),
  list(
    name = "endpoint-first", dash = "13", width_pt = 1.6,
    formula = quote(phantom() %~~% paste("(", n - 1, ")", K^2)),
    fun = ops_endpoint
  )
)

scientific_label <- function(value) {
  exponent <- floor(log10(value))
  mantissa <- value / 10^exponent
  mantissa_text <- sprintf("%.1f", mantissa)
  expr <- bquote(paste(.(mantissa_text), " " %*% " ", 10^.(exponent)))
  expr
}

plot_x <- function(n) {
  span <- plot_right - plot_left
  value <- plot_left + (n - n_min) / (n_max - n_min) * span
  value
}

plot_y <- function(value) {
  fraction <- (log10(value) - log_min) / (log_max - log_min)
  height <- plot_bottom - plot_top
  result <- plot_bottom - fraction * height
  result
}

draw_count_table <- function() {
  draw_prose(40, 760, "order", 11.5, col_ink)
  draw_prose(186, 760, "operations", 11.5, col_ink)
  draw_prose(330, 760, "n = 5, K = 10", 11.5, col_ink)
  y_row <- 786
  for (entry in count_rows) {
    line_width <- lwd_from_pt(entry$width_pt)
    segments(40, y_row - 4, 74, y_row - 4,
      col = col_ink, lwd = line_width,
      lty = entry$dash
    )
    draw_prose(84, y_row, entry$name, 12.5, col_caption)
    draw_math(186, y_row, entry$formula, 14, col_ink, 0)
    worked_value <- entry$fun(5, k_states)
    worked_label <- scientific_label(worked_value)
    draw_math(330, y_row, worked_label, 12.5, col_caption, 0, "sans")
    y_row <- y_row + 26
  }
  invisible(NULL)
}

draw_count_plot <- function() {
  axis_width <- lwd_from_pt(1.2)
  segments(plot_left, plot_top, plot_left, plot_bottom,
    col = col_ink,
    lwd = axis_width
  )
  segments(plot_left, plot_bottom, plot_right, plot_bottom,
    col = col_ink,
    lwd = axis_width
  )
  tick_width <- lwd_from_pt(1)
  for (exponent in c(2, 4, 6, 8, 10)) {
    y_tick <- plot_y(10^exponent)
    segments(plot_left - 5, y_tick, plot_left, y_tick,
      col = col_ink,
      lwd = tick_width
    )
    tick_label <- bquote(10^.(exponent))
    draw_math(
      plot_left - 9, y_tick + 4, tick_label, 11,
      col_ink, 1, "sans"
    )
  }
  for (n_value in c(4, 6, 8, 10)) {
    x_tick <- plot_x(n_value)
    segments(x_tick, plot_bottom, x_tick, plot_bottom + 5,
      col = col_ink,
      lwd = tick_width
    )
    draw_prose(
      x_tick, plot_bottom + 18, as.character(n_value), 11,
      col_ink, 0.5
    )
  }
  draw_prose(plot_left - 58, (plot_top + plot_bottom) / 2, "operations",
    11.5, col_ink,
    align = 0.5, rotation = 90
  )
  x_axis_centre <- (plot_left + plot_right) / 2
  draw_math(
    x_axis_centre, plot_bottom + 34, quote(n), 13, col_ink,
    0.5
  )
  x_worked <- plot_x(5)
  segments(x_worked, plot_top, x_worked, plot_bottom,
    col = col_rule,
    lwd = tick_width
  )
  n_grid <- n_min:n_max
  for (entry in count_rows) {
    y_values <- entry$fun(n_grid, k_states)
    x_curve <- plot_x(n_grid)
    y_curve <- plot_y(y_values)
    line_width <- lwd_from_pt(entry$width_pt)
    lines(x_curve, y_curve,
      col = col_ink, lwd = line_width,
      lty = entry$dash
    )
    y_worked <- plot_y(entry$fun(5, k_states))
    point_size <- cex_from_pt(8)
    points(x_worked, y_worked, pch = 16, col = col_ink, cex = point_size)
  }
  draw_prose(
    x_worked + 6, plot_top + 14, "worked example", 10.5,
    col_gone_text
  )
  y_brute <- plot_y(ops_brute(n_max - 1, k_states)) + 26
  draw_prose(plot_right - 6, y_brute, "brute force", 11, col_caption, 1)
  y_interior <- plot_y(ops_interior(n_max, k_states)) - 9
  draw_prose(
    plot_right - 6, y_interior, "interior-first", 11, col_caption,
    1
  )
  y_endpoint <- plot_y(ops_endpoint(n_max, k_states)) + 16
  draw_prose(
    plot_right - 6, y_endpoint, "endpoint-first", 11, col_caption,
    1
  )
  invisible(NULL)
}

# ---- assemble -----------------------------------------------------------

draw_figure <- function() {
  par(mar = rep(0, 4), ps = 12, xpd = NA)
  plot.new()
  plot.window(
    xlim = c(0, fig_width_px),
    ylim = c(fig_height_px, 0),
    xaxs = "i",
    yaxs = "i",
    asp = 1
  )

  head_a_1 <- paste0(
    "(a)  The marginal likelihood sums a product of factors over every ",
    "joint configuration of the latent variables. Two"
  )
  head_a_2 <- paste0(
    "variables are joined in the variable graph whenever they appear in ",
    "the same factor."
  )
  draw_prose(40, 26, head_a_1, 13.5, col_ink)
  draw_prose(40, 44, head_a_2, 13.5, col_ink)

  setup_lhs <- quote(paste(L, "  =  "))
  setup_index <- quote(paste(x[1], ",", ldots, ",", x[5]))
  setup_rhs_left <- join_two(
    pair_factor("f", 12, 1, 2),
    pair_factor("f", 23, 2, 3)
  )
  setup_rhs_right <- join_two(
    pair_factor("f", 34, 3, 4),
    pair_factor("f", 45, 4, 5)
  )
  setup_rhs <- join_two(setup_rhs_left, setup_rhs_right)
  draw_equation(560, y_setup_eq, setup_lhs, setup_index, setup_rhs, 16)

  shift <- 560 - node_x(1, 3)
  panel_left[1] <<- panel_left[1] + shift
  draw_graph(1, frame_setup, y_setup_node)
  panel_left[1] <<- panel_left[1] - shift

  setup_note <- bquote(
    paste(
      "chain length ", n, " = 5; each ", x[i], " takes ", K,
      " = 10 possible values; brute force enumerates all ", K^n,
      " configurations: "
    ) %~~%
      paste(4 * K^5, " operations")
  )
  draw_math(
    560, y_setup_note, setup_note, 12.5, col_caption, 0.5,
    "sans"
  )

  draw_rule(200)

  head_b_1 <- bquote(paste(
    "(b)  Interior-first order (", x[2], ", ", x[4], ", ", x[3],
    "). Each sum couples two surviving neighbours, so its result is a ",
    "new pairwise factor: a fill-in"
  ))
  head_b_2 <- bquote(paste(
    "edge. Maximum clique 3, total ", 3 * K^3, " + ", K^2, "."
  ))
  draw_math(40, 224, head_b_1, 13.5, col_ink, 0, "sans")
  draw_math(40, 248, head_b_2, 13.5, col_ink, 0, "sans")
  draw_row(row_b, y_b_eq, y_b_node, y_b_caption, accent = col_purple)

  draw_rule(454)

  head_c_1 <- bquote(paste(
    "(c)  Endpoint-first order (", x[1], ", ", x[2], ", ", x[3],
    "). Each sum touches a single surviving neighbour, so its result is ",
    "univariate and no edge is"
  ))
  head_c_2 <- bquote(paste(
    "added. Maximum clique 2, total ", 4 * K^2, " + K, which is the ",
    "forward recursion."
  ))
  draw_math(40, 478, head_c_1, 13.5, col_ink, 0, "sans")
  draw_math(40, 502, head_c_2, 13.5, col_ink, 0, "sans")
  draw_row(row_c, y_c_eq, y_c_node, y_c_caption, accent = col_purple)

  draw_rule(708)

  head_d <- bquote(paste(
    "(d)  Operation counts as chain length ", n,
    " grows (K states each; panels a-c use ", n, " = 5)."
  ))
  draw_math(40, 732, head_d, 13.5, col_ink, 0, "sans")
  draw_count_table()
  draw_count_plot()

  foot_1 <- paste0(
    "The two elimination orders differ by a factor of K; both differ from"
  )
  foot_2 <- "enumeration by a margin that grows without bound in n."
  draw_prose(40, 892, foot_1, 11.5, col_ink)
  draw_prose(40, 908, foot_2, 11.5, col_ink)
  invisible(NULL)
}

# ---- device -------------------------------------------------------------

# Scaling width/height and res together keeps the physical page size fixed
# (inches = pixels / res is unchanged) while raising the pixel density, so
# every hand-placed coordinate above stays correct at any render_scale.
render_scale <- 3

out_path <- "figures/sequential_reduction.png"
png(out_path,
  width = fig_width_px * render_scale,
  height = fig_height_px * render_scale,
  pointsize = 12,
  bg = "white",
  res = 96 * render_scale
)
draw_figure()
dev.off()
