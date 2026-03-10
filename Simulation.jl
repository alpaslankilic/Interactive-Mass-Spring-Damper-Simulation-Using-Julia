using DifferentialEquations
import GLMakie

# ── ODE ───────────────────────────────────────────────────────────────────────
function smd!(du, u, p, t)
    x, v = u
    m, b, k = p
    du[1] = v
    du[2] = -(b*v + k*x) / m
end

# ── Spring geometry ───────────────────────────────────────────────────────────
function make_vspring(y_top, y_bot, x_ctr; n_half=14, amp=0.055)
    lead   = (y_top - y_bot) * 0.08
    coil_h = (y_top - y_bot - 2*lead) / n_half
    xs = Float64[x_ctr, x_ctr]
    ys = Float64[y_top, y_top - lead]
    for i in 1:n_half
        push!(xs, x_ctr + (isodd(i) ? amp : -amp))
        push!(ys, y_top - lead - i * coil_h)
    end
    push!(xs, x_ctr); push!(ys, y_bot)
    return xs, ys
end

# ── Constants ─────────────────────────────────────────────────────────────────
const WIN     = 6.0
const SPR_X   = -0.14
const DMP_X   =  0.14
const CYL_HW  =  0.050
const STEP_DT =  0.02

# ── Figure ────────────────────────────────────────────────────────────────────
fig = GLMakie.Figure(size=(1160, 740), backgroundcolor=:grey93)

# ── Slider grid ───────────────────────────────────────────────────────────────
sg = GLMakie.SliderGrid(fig[2, 1:2],
    (label = "m  (kg)",    range = 0.1:0.1:5.0,  startvalue = 1.0),
    (label = "b  (N·s/m)", range = 0.0:0.1:10.0,  startvalue = 0.5),
    (label = "k  (N/m)",   range = 0.5:0.5:10.0, startvalue = 5.0),
    (label = "x0 (m)",     range = -1.0:0.05:1.0, startvalue = 0.4);
    tellwidth = false)

m_s  = sg.sliders[1].value
b_s  = sg.sliders[2].value
k_s  = sg.sliders[3].value
x0_s = sg.sliders[4].value

# ── Button row ────────────────────────────────────────────────────────────────
btn_row = fig[3, 1:2] = GLMakie.GridLayout(tellwidth = false)
inner   = btn_row[1, 1] = GLMakie.GridLayout(tellwidth = false)
GLMakie.colsize!(btn_row, 1, GLMakie.Relative(1.0))

btn_play    = GLMakie.Button(inner[1, 1];
    label = " > ", fontsize = 22, width = 64, height = 46,
    buttoncolor = :grey85, buttoncolor_hover = :lightgreen,
    buttoncolor_active = :green, strokewidth = 1)

btn_pause   = GLMakie.Button(inner[1, 2];
    label = " || ", fontsize = 22, width = 64, height = 46,
    buttoncolor = :grey85, buttoncolor_hover = :lightyellow,
    buttoncolor_active = :yellow, strokewidth = 1)

btn_restart = GLMakie.Button(inner[1, 3];
    label = " |< ", fontsize = 22, width = 64, height = 46,
    buttoncolor = :grey85, buttoncolor_hover = :lightblue,
    buttoncolor_active = :dodgerblue, strokewidth = 1)

btn_close   = GLMakie.Button(inner[1, 4];
    label = " X ", fontsize = 22, width = 64, height = 46,
    buttoncolor = :grey85,
    buttoncolor_hover   = GLMakie.RGBf(1.0, 0.55, 0.55),
    buttoncolor_active  = GLMakie.RGBf(0.85, 0.1, 0.1),
    strokewidth = 1)

GLMakie.colgap!(inner, 10)

# ── Simulation state ──────────────────────────────────────────────────────────
is_running = Ref(true)
do_close   = Ref(false)

# ── Integrator factory ────────────────────────────────────────────────────────
function make_integ(m, b, k, x0)
    prob = ODEProblem(smd!, [x0, 0.0], (0.0, 1e10), (m, b, k))
    init(prob, Tsit5(); reltol=1e-8, abstol=1e-10,
         dtmax=STEP_DT, save_everystep=false)
end

integ = Ref(make_integ(m_s[], b_s[], k_s[], x0_s[]))

# ── Rolling history buffers ───────────────────────────────────────────────────
t_buf = Float64[0.0]
x_buf = Float64[x0_s[]]

# ── Live observables ──────────────────────────────────────────────────────────
x_now  = GLMakie.Observable(x0_s[])
hist_t = GLMakie.Observable(copy(t_buf))
hist_x = GLMakie.Observable(copy(x_buf))
cur_pt = GLMakie.Observable([GLMakie.Point2f(0.0, x0_s[])])

# y_max driven by abs(x0) so negative values scale correctly
y_max_obs = GLMakie.@lift max(abs($(x0_s)) * 1.35, 0.15)
CEIL_Y    = GLMakie.@lift  $(y_max_obs) * 1.08
BLK_HH    = GLMakie.@lift  $(y_max_obs) * 0.062
CYL_H     = GLMakie.@lift  $(y_max_obs) * 0.28
cyl_bot   = GLMakie.@lift  $(CEIL_Y) - $(CYL_H)
blk_top   = GLMakie.@lift  $x_now + $(BLK_HH)
blk_bot   = GLMakie.@lift  $x_now - $(BLK_HH)

# ── Reset / restart ───────────────────────────────────────────────────────────
function restart_sim!()
    integ[] = make_integ(m_s[], b_s[], k_s[], x0_s[])
    empty!(t_buf); empty!(x_buf)
    push!(t_buf, 0.0); push!(x_buf, x0_s[])
    x_now[]  = x0_s[]
    hist_t[] = copy(t_buf)
    hist_x[] = copy(x_buf)
    cur_pt[] = [GLMakie.Point2f(0.0, x0_s[])]
    is_running[] = true
end

# ── Button callbacks ──────────────────────────────────────────────────────────
GLMakie.on(btn_play.clicks)    do _; is_running[] = true        end
GLMakie.on(btn_pause.clicks)   do _; is_running[] = false       end
GLMakie.on(btn_restart.clicks) do _; restart_sim!()             end
GLMakie.on(btn_close.clicks)   do _
    is_running[] = false
    do_close[]   = true
    GLMakie.closeall()
end

# Restart when any slider moves
for s in [m_s, b_s, k_s, x0_s]
    GLMakie.on(s) do _; restart_sim!() end
end

# ── Axes ──────────────────────────────────────────────────────────────────────
ax_g = GLMakie.Axis(fig[1, 1];
    xlabel            = "Time  (s)",
    ylabel            = "Displacement  (m)",
    backgroundcolor   = :white,
    xgridvisible      = true,
    ygridvisible      = true,
    rightspinevisible = false,
    topspinevisible   = false)

ax_m = GLMakie.Axis(fig[1, 2];
    backgroundcolor    = :white,
    xgridvisible       = false,
    ygridvisible       = false,
    yticksvisible      = false,
    yticklabelsvisible = false,
    xticksvisible      = false,
    xticklabelsvisible = false,
    leftspinevisible   = false,
    rightspinevisible  = false,
    topspinevisible    = false,
    bottomspinevisible = false)

GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.67))
GLMakie.colsize!(fig.layout, 2, GLMakie.Relative(0.33))

GLMakie.on(y_max_obs) do ym
    GLMakie.ylims!(ax_g, -ym,      ym)
    GLMakie.ylims!(ax_m, -ym*1.25, ym*1.30)
    GLMakie.xlims!(ax_m, -0.38, 0.38)
    GLMakie.xlims!(ax_g,  0.0,  WIN)
end
GLMakie.notify(y_max_obs)

# ── Graph elements ────────────────────────────────────────────────────────────
GLMakie.hlines!(ax_g, [0.0]; color=(:grey55, 0.7), linewidth=1.5, linestyle=:dash)
GLMakie.lines!(ax_g, hist_t, hist_x; color=:steelblue, linewidth=2.8)
GLMakie.scatter!(ax_g, cur_pt;
    color=:crimson, markersize=13, strokewidth=1.5, strokecolor=:white)

# ── Mechanical panel ──────────────────────────────────────────────────────────
# Ceiling bar
GLMakie.lines!(ax_m, [-0.36, 0.36],
    GLMakie.@lift([$(CEIL_Y), $(CEIL_Y)]); color=:black, linewidth=3.5)

# Ceiling hatch marks
GLMakie.linesegments!(ax_m,
    GLMakie.@lift(begin
        pts = GLMakie.Point2f[]
        cy  = $(CEIL_Y); off = $(y_max_obs) * 0.07
        for xh in -0.34:0.09:0.34
            push!(pts, GLMakie.Point2f(xh,         cy))
            push!(pts, GLMakie.Point2f(xh + 0.055, cy + off))
        end
        pts
    end); color=:grey40, linewidth=1.3)

# Equilibrium dashed line in mechanical panel
GLMakie.hlines!(ax_m, [0.0]; color=(:grey55, 0.5), linewidth=1.0, linestyle=:dash)
# Initial displacement line (x₀)
GLMakie.hlines!(ax_m,
    GLMakie.@lift([$(x0_s)]);
    color = (:crimson, 0.7),
    linewidth = 1.5,
    linestyle = :dashdot)
GLMakie.text!(ax_m,
    GLMakie.@lift(GLMakie.Point2f(0.26, $(x0_s)));
    text = GLMakie.@lift("x₀ = $(round($(x0_s), digits=2))"),
    align = (:left, :top),
    fontsize = 16,
    color = :crimson)
# Spring
GLMakie.lines!(ax_m,
    GLMakie.@lift(make_vspring($(CEIL_Y), $(blk_top), SPR_X)[1]),
    GLMakie.@lift(make_vspring($(CEIL_Y), $(blk_top), SPR_X)[2]);
    color=:royalblue, linewidth=2.6)

# ── Improved damper ───────────────────────────────────────────────────────────

# piston position inside the cylinder
# it moves with the mass, but is clamped so it never exits the cylinder
piston_y = GLMakie.@lift begin
    margin = $(CYL_H) * 0.15

    # map mass displacement into the cylinder interior
    raw = $(cyl_bot) + 0.20 * $(CYL_H) + (($(blk_top) - 0.0) / $(y_max_obs)) * (0.60 * $(CYL_H))

    clamp(raw, $(cyl_bot) + margin, $(CEIL_Y) - margin)
end

# small rod segment inside the cylinder (from bottom opening to piston)
rod_in_top = piston_y
rod_in_bot = GLMakie.@lift $(cyl_bot) + $(CYL_H) * 0.06

# visible rod outside the cylinder (from mass to cylinder bottom)
GLMakie.lines!(ax_m, [DMP_X, DMP_X],
    GLMakie.@lift([$(blk_top), $(cyl_bot)]);
    color=:grey25, linewidth=2.8)

# cylinder body
GLMakie.poly!(ax_m,
    GLMakie.@lift(GLMakie.Point2f[
        (DMP_X - CYL_HW, $(cyl_bot)),
        (DMP_X + CYL_HW, $(cyl_bot)),
        (DMP_X + CYL_HW, $(CEIL_Y)),
        (DMP_X - CYL_HW, $(CEIL_Y))
    ]);
    color=(:lightgrey, 0.88), strokecolor=:grey35, strokewidth=2.0)

# internal rod inside the cylinder
GLMakie.lines!(ax_m, [DMP_X, DMP_X],
    GLMakie.@lift([$(rod_in_bot), $(rod_in_top)]);
    color=:grey25, linewidth=2.4)

# piston head (moving)
GLMakie.linesegments!(ax_m,
    GLMakie.@lift([
        GLMakie.Point2f(DMP_X - CYL_HW + 0.008, $(piston_y)),
        GLMakie.Point2f(DMP_X + CYL_HW - 0.008, $(piston_y))
    ]);
    color=:grey30, linewidth=4.6)

# bottom seal / guide where the rod enters the cylinder
GLMakie.linesegments!(ax_m,
    GLMakie.@lift([
        GLMakie.Point2f(DMP_X - CYL_HW * 0.24, $(cyl_bot)),
        GLMakie.Point2f(DMP_X + CYL_HW * 0.24, $(cyl_bot))
    ]);
    color=:grey20, linewidth=3.2)

# k / b labels WITH VALUES
lbl_y = GLMakie.@lift $(CEIL_Y) - $(y_max_obs) * 0.17

GLMakie.text!(ax_m,
    GLMakie.@lift(GLMakie.Point2f(SPR_X -0.19, $(lbl_y)));
    text = GLMakie.@lift("k = $(round($(k_s), digits=2))"),
    color = :royalblue,
    fontsize = 17,
    font = :bold)

GLMakie.text!(ax_m,
    GLMakie.@lift(GLMakie.Point2f(DMP_X + 0.07, $(lbl_y)));
    text = GLMakie.@lift("b = $(round($(b_s), digits=2))"),
    color = :grey30,
    fontsize = 17,
    font = :bold)

# Mass block
GLMakie.poly!(ax_m,
    GLMakie.@lift(GLMakie.Point2f[
        (-0.24, $(blk_bot)), ( 0.24, $(blk_bot)),
        ( 0.24, $(blk_top)), (-0.24, $(blk_top))]);
    color=(:steelblue, 0.90), strokecolor=:navy, strokewidth=2.5)

GLMakie.text!(ax_m,
    GLMakie.@lift(GLMakie.Point2f(0.0, ($(blk_top) + $(blk_bot)) * 0.5));
    text = GLMakie.@lift("m = $(round($(m_s), digits=2))"),
    align = (:center, :center),
    color = :white,
    fontsize = 20,
    font = :bold)

# ── Display & animation loop ──────────────────────────────────────────────────
display(fig)

fps_a    = 30
frame_dt = 1.0 / fps_a

while GLMakie.isopen(fig.scene) && !do_close[]
    if is_running[]
        step!(integ[], STEP_DT, true)
        t_i = integ[].t
        x_i = integ[].u[1]

        push!(t_buf, t_i)
        push!(x_buf, x_i)

        # Trim buffer to last WIN*1.5 s
        while length(t_buf) > 1 && (t_buf[end] - t_buf[1]) > WIN * 1.5
            popfirst!(t_buf); popfirst!(x_buf)
        end

        x_now[]  = x_i
        hist_t[] = t_buf
        hist_x[] = x_buf
        cur_pt[] = [GLMakie.Point2f(t_i, x_i)]

        t_lo = max(t_buf[1], t_i - WIN)
        GLMakie.xlims!(ax_g, t_lo, t_lo + WIN)
    end

    sleep(frame_dt)
end
