import ast, math
import ipywidgets as W
from IPython.display import display, clear_output
from sage.all import Link, Graph, QQ  # assumes Sage kernel


# -------------------------
# UI helper widgets
# -------------------------

def helpbox(md):
    return W.HTML(
        f"""
        <div style="padding:10px; border:1px solid #444; border-radius:8px;
                    background:#1f1f1f; color:#ddd; line-height:1.3;">
        {md}
        </div>
        """
    )


# -------------------------
# Plotting / PD helpers
# -------------------------

def pd_for_sage_plot(pd):
    if not pd:
        return [], {}
    labels = sorted(set(abs(x) for cr in pd for x in cr))
    mp = {lab: i + 1 for i, lab in enumerate(labels)}
    pd_plot = [[mp[abs(x)] for x in cr] for cr in pd]
    return pd_plot, mp


def _all_strands(link):
    pd = link.pd_code()
    if not pd:
        return [-1]
    return sorted(set(x for cr in pd for x in cr))


def _loop_strands(link):
    pd = link.pd_code()
    if not pd:
        return []
    loops = []
    strands = _all_strands(link)
    for s in strands:
        for cr in pd:
            if sum(1 for x in cr if x == s) >= 2:
                loops.append(s)
                break
    return sorted(set(loops))


def pd_connected_components(pd):
    """
    Components of the PD by connectivity through shared strand labels.
    Returns list[set[int]] of crossing indices.
    """
    n = len(pd)
    if n == 0:
        return []

    occ = {}
    for i, cr in enumerate(pd):
        for s in cr:
            occ.setdefault(s, set()).add(i)

    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for s, idxs in occ.items():
        idxs = list(idxs)
        if len(idxs) >= 2:
            base = idxs[0]
            for j in idxs[1:]:
                union(base, j)

    comps = {}
    for i in range(n):
        comps.setdefault(find(i), set()).add(i)

    out = list(comps.values())
    out.sort(key=lambda C: (-len(C), min(C)))
    return out


def strand_graph_plot_stable(link, title=None):
    """
    Strand graph with deterministic layout:
    each PD-connected component placed on a separate circle.
    Vertices = crossing indices, edge labels = strand labels.
    """
    pd = link.pd_code()
    if not pd:
        G = Graph(loops=True, multiedges=True)
        return G.plot(title=title or "Unknot (no crossings)")

    occ = {}
    for i, cr in enumerate(pd):
        for s in cr:
            occ.setdefault(s, []).append(i)

    G = Graph(loops=True, multiedges=True)
    G.add_vertices(list(range(len(pd))))

    for s, idxs in occ.items():
        idxs = list(set(idxs))
        if len(idxs) == 2:
            u, v = idxs
            G.add_edge(u, v, str(s))
        elif len(idxs) == 1:
            u = idxs[0]
            G.add_edge(u, u, str(s))  # loop edge
        else:
            base = idxs[0]
            for v in idxs[1:]:
                G.add_edge(base, v, str(s))

    comps = pd_connected_components(pd)

    pos = {}
    x_offset = 0.0
    for comp in comps:
        comp = sorted(comp)
        m = len(comp)
        r = 1.2 if m > 1 else 0.15
        for k, v in enumerate(comp):
            ang = 2 * math.pi * k / max(1, m)
            pos[v] = (x_offset + r * math.cos(ang), r * math.sin(ang))
        x_offset += 3.5

    return G.plot(
        pos=pos,
        vertex_labels=True,
        edge_labels=True,
        title=title or "Strand graph (stable layout: edges=strand labels; vertices=crossing indices)"
    )


# -------------------------
# Main UI
# -------------------------

def launch_movie_ui(Movie, start_pd=None, starting_qdeg=1, ring=QQ):
    if start_pd is None:
        start_pd = []
    movie = Movie(Link(start_pd), starting_qdeg=starting_qdeg, ring=ring)

    out_main = W.Output()
    out_side = W.Output()
    out_log  = W.Output()

    play = W.Play(interval=350, value=0, min=0, max=0, step=1, description="Play")
    frame = W.IntSlider(value=0, min=0, max=0, step=1, description="Frame", continuous_update=False)
    W.jslink((play, 'value'), (frame, 'value'))

    show_strand_graph = W.Checkbox(value=True, description="Show strand graph (numbered)")
    show_pd = W.Checkbox(value=True, description="Show PD + strand list")
    plot_components_separately = W.Checkbox(value=True, description="Plot components separately")

    pd_text = W.Textarea(
        value=str(start_pd),
        description="Start PD",
        layout=W.Layout(width="520px", height="90px")
    )
    btn_load = W.Button(description="Load / Reset", button_style="warning")
    btn_jump_last = W.Button(description="Jump to last frame")
    btn_clear_log = W.Button(description="Clear log")

    # ---- widgets (renamed) ----
    dd_twist_strand = W.Dropdown(description="Target strand", options=[-1], value=-1)
    tg_twist_orient = W.ToggleButtons(description="Loop orient", options=[1, -1], value=1)
    tg_twist_type   = W.ToggleButtons(description="Over/under", options=[1, -1], value=1)
    btn_twist = W.Button(description="Apply R1 (twist)", button_style="primary")

    dd_untwist_loop = W.Dropdown(description="Loop strand", options=[], value=None)
    btn_untwist = W.Button(description="Apply R1⁻¹ (untwist)", button_style="primary")

    dd_poke_s1 = W.Dropdown(description="Moving strand", options=[-1], value=-1)
    dd_poke_s2 = W.Dropdown(description="Stationary strand", options=[-1], value=-1)
    tg_poke_parallel = W.ToggleButtons(description="Parallel?", options=[-1, 1], value=-1)
    tg_poke_over     = W.ToggleButtons(description="Moving over?", options=[1, -1], value=1)
    btn_poke = W.Button(description="Apply R2 (poke)", button_style="primary")

    dd_unpoke_s1 = W.Dropdown(description="Moving strand", options=[-1], value=-1)
    dd_unpoke_s2 = W.Dropdown(description="Stationary strand", options=[-1], value=-1)
    btn_unpoke = W.Button(description="Apply R2⁻¹ (unpoke)", button_style="primary")

    dd_slide_lms = W.Dropdown(description="Left moving", options=[-1], value=-1)
    tg_slide_om  = W.ToggleButtons(description="Moving dir", options=[1, -1], value=1)
    tg_slide_ol  = W.ToggleButtons(description="Left dir", options=[1, -1], value=1)
    tg_slide_or  = W.ToggleButtons(description="Right dir", options=[1, -1], value=1)
    btn_slide = W.Button(description="Apply R3 (slide)", button_style="primary")

    dd_saddle_s1 = W.Dropdown(description="Strand A", options=[-1], value=-1)
    dd_saddle_s2 = W.Dropdown(description="Strand B", options=[-1], value=-1)
    btn_saddle = W.Button(description="Apply saddle", button_style="primary")

    btn_birth = W.Button(description="Birth component", button_style="primary")
    dd_death_loop = W.Dropdown(description="Loop strand", options=[], value=None)
    btn_death = W.Button(description="Death component", button_style="primary")

    perm_text = W.Textarea(
        value="{ }",
        description="Permutation",
        layout=W.Layout(width="520px", height="70px")
    )
    btn_perm = W.Button(description="Apply permutation", button_style="primary")

    # ---- tabs with explanations ----
    tab = W.Tab(children=[
        W.VBox([
            helpbox("""
            <b>R1: Add a twist</b><br>
            Adds one crossing that creates a small loop/kink on the chosen strand.<br><br>
            <b>Target strand</b>: pick an existing strand label.<br>
            <b>Loop orient</b> (±1): chooses convention for whether the 0-resolution contains the loop.<br>
            <b>Over/under</b> (±1): +1 = crossing treated as overstrand, −1 = understrand.
            """),
            W.HBox([dd_twist_strand, tg_twist_orient, tg_twist_type]),
            btn_twist
        ]),
        W.VBox([
            helpbox("""
            <b>R1: Remove a twist</b><br>
            Removes a loop-crossing (inverse of R1).<br><br>
            <b>Loop strand</b>: must be a label that appears twice in a single crossing (dropdown only shows candidates).
            """),
            dd_untwist_loop,
            btn_untwist
        ]),
        W.VBox([
            helpbox("""
            <b>R2: Poke (add two crossings)</b><br>
            Creates a Reidemeister-II pair by pushing one strand around another.<br><br>
            <b>Moving strand</b>: strand that does the poking.<br>
            <b>Stationary strand</b>: strand being poked around.<br>
            <b>Parallel?</b>: +1 if oriented strands are parallel, −1 otherwise.<br>
            <b>Moving over?</b>: +1 if moving strand goes over, −1 if under.
            """),
            W.HBox([dd_poke_s1, dd_poke_s2]),
            W.HBox([tg_poke_parallel, tg_poke_over]),
            btn_poke
        ]),
        W.VBox([
            helpbox("""
            <b>R2: Unpoke (remove two crossings)</b><br>
            Inverse of <i>Poke</i>. Use when you have a clean R2 pair to remove.<br><br>
            <b>Inputs</b>: pick the same conceptual roles as in <i>Poke</i>:
            the moving strand (the one that loops around) and the stationary strand.
            """),
            W.HBox([dd_unpoke_s1, dd_unpoke_s2]),
            btn_unpoke
        ]),
        W.VBox([
            helpbox("""
            <b>R3: Slide</b><br>
            Performs a Reidemeister-III move (strand slides over/under a crossing).<br><br>
            <b>Left moving</b>: pick the left segment of the moving strand in top-to-bottom convention.<br>
            <b>Moving/Left/Right dir</b> (±1): orientations as in <code>slide</code> docstring.
            """),
            W.HBox([dd_slide_lms, tg_slide_om]),
            W.HBox([tg_slide_ol, tg_slide_or]),
            btn_slide
        ]),
        W.VBox([
            helpbox("""
            <b>Saddle</b><br>
            Performs a saddle cobordism between two strands (merge/split depending on local picture).<br><br>
            <b>Strand A/B</b>: strands where you attach the band (must exist in the PD).<br>
            <b>Note</b>: Code decrements <code>last_degree</code> by 1.
            """),
            W.HBox([dd_saddle_s1, dd_saddle_s2]),
            btn_saddle
        ]),
        W.VBox([
            helpbox("""
            <b>Birth / Death</b><br>
            Adds/removes a separated component using “twisted unknot” convention.<br><br>
            <b>Birth</b>: adds a new twisted unknot component (may first twist if starting from unknot).<br>
            <b>Death</b>: removes a twisted unknot component; choose a loop-strand on that component.<br>
            <b>Tip</b>: for sign-sensitive work, follows warning about which loop label to kill.
            """),
            btn_birth,
            W.HTML("<hr>"),
            dd_death_loop,
            btn_death
        ]),
        W.VBox([
            helpbox("""
            <b>Permute crossings</b><br>
            Reorders crossing indices and composes the induced chain isomorphism into the last map.<br><br>
            Enter a dict <code>{i:j, ...}</code> meaning “send crossing index i → j”.<br>
            Example: <code>{0:2, 2:0}</code>.
            """),
            perm_text,
            btn_perm
        ]),
    ])

    for i, name in enumerate(["R1: Twist", "R1: Untwist", "R2: Poke", "R2: Unpoke", "R3: Slide", "Saddle", "Birth/Death", "Permute"]):
        tab.set_title(i, name)

    # -------------------------
    # Internal UI logic
    # -------------------------

    def refresh_dropdowns():
        L = movie.links[-1]
        strands = _all_strands(L)
        loops = _loop_strands(L)

        for dd in [dd_twist_strand, dd_poke_s1, dd_poke_s2, dd_unpoke_s1, dd_unpoke_s2,
                   dd_slide_lms, dd_saddle_s1, dd_saddle_s2]:
            dd.options = strands
            if dd.value not in strands:
                dd.value = strands[0]

        dd_untwist_loop.options = loops
        dd_death_loop.options = loops
        if loops:
            if dd_untwist_loop.value not in loops:
                dd_untwist_loop.value = loops[0]
            if dd_death_loop.value not in loops:
                dd_death_loop.value = loops[0]
        else:
            dd_untwist_loop.value = None
            dd_death_loop.value = None

    def render():
        n = len(movie.links)
        play.max = max(0, n - 1)
        frame.max = max(0, n - 1)

        idx = int(frame.value)
        L = movie.links[idx]
        pd0 = L.pd_code()

        with out_main:
            clear_output(wait=True)
            display(W.HTML(
                f"<b>Frames:</b> {n} &nbsp;&nbsp; "
                f"<b>Current frame:</b> {idx} &nbsp;&nbsp; "
                f"<b>Current q-degree:</b> {movie.last_degree} &nbsp;&nbsp; "
                f"<b>Min strand label:</b> {movie.current_min_index}"
            ))

            # Component info (helps interpret plots)
            comps = pd_connected_components(pd0) if pd0 else []
            display(W.HTML(f"<b>PD components:</b> {len(comps)}"))
            if comps:
                for ci, C in enumerate(comps):
                    display(W.HTML(f"<pre>Component {ci}: crossings {sorted(C)}</pre>"))

            # Plot relabel mapping for clarity
            if pd0:
                _, mp = pd_for_sage_plot(pd0)
                display(W.HTML(f"<b>Plot relabel mapping</b> (abs(label) → plotlabel):<pre>{mp}</pre>"))

            # Main plot (robust)
            try:
                if plot_components_separately.value and pd0:
                    for ci, C in enumerate(comps):
                        sub_pd = [pd0[i] for i in sorted(C)]
                        pd_plot_sub, _ = pd_for_sage_plot(sub_pd)
                        display(Link(pd_plot_sub).plot(title=f"Component {ci} (frame {idx}, plot labels renumbered)"))
                else:
                    pd_plot, _ = pd_for_sage_plot(pd0)
                    display(Link(pd_plot).plot(title=f"Link diagram (frame {idx})  (plot labels renumbered)"))
            except Exception as e:
                display(W.HTML(f"<pre>Sage Link.plot() failed: {e}\nShowing strand-graph fallback instead.</pre>"))
                try:
                    display(strand_graph_plot_stable(L, title="Fallback strand graph (stable layout)"))
                except Exception as e2:
                    display(W.HTML(f"<pre>Fallback strand graph also failed: {e2}</pre>"))

        with out_side:
            clear_output(wait=True)
            if show_strand_graph.value:
                try:
                    display(strand_graph_plot_stable(L))
                except Exception as e:
                    display(W.HTML(f"<pre>Strand-graph plot error: {e}</pre>"))

            if show_pd.value:
                display(W.HTML("<b>PD code</b>:"))
                display(W.HTML(f"<pre>{pd0}</pre>"))
                display(W.HTML("<b>Strands present</b>:"))
                display(W.HTML(f"<pre>{_all_strands(L)}</pre>"))
                if pd0:
                    display(W.HTML("<b>Loop strands (candidates for untwist/death)</b>:"))
                    display(W.HTML(f"<pre>{_loop_strands(L)}</pre>"))

    def log(msg):
        with out_log:
            print(msg)

    def jump_to_last(_=None):
        frame.value = max(0, len(movie.links) - 1)
        render()

    def on_load(_):
        nonlocal movie
        try:
            new_pd = ast.literal_eval(pd_text.value.strip())
            if not isinstance(new_pd, list):
                raise ValueError("Start PD must be a Python list, e.g. [] or [[-1,-2,-2,-1]]")
            movie = Movie(Link(new_pd), starting_qdeg=starting_qdeg, ring=ring)
            frame.value = 0
            refresh_dropdowns()
            render()
            log("Loaded/reset movie.")
        except Exception as e:
            log(f"[load error] {e}")

    btn_load.on_click(on_load)
    btn_jump_last.on_click(jump_to_last)
    btn_clear_log.on_click(lambda _: out_log.clear_output())

    def apply_move(fn_name, *args):
        try:
            getattr(movie, fn_name)(*args, print_pd=False)
            refresh_dropdowns()
            jump_to_last()
            log(f"Applied {fn_name}{args}")
        except Exception as e:
            log(f"[{fn_name} error] {e}")

    btn_twist.on_click(lambda _:
        apply_move("twist", dd_twist_strand.value, tg_twist_orient.value, tg_twist_type.value)
    )
    btn_untwist.on_click(lambda _:
        apply_move("untwist", dd_untwist_loop.value) if dd_untwist_loop.value is not None else log("[untwist] no loop strands available")
    )
    btn_poke.on_click(lambda _:
        apply_move("poke", dd_poke_s1.value, dd_poke_s2.value, tg_poke_parallel.value, tg_poke_over.value)
    )
    btn_unpoke.on_click(lambda _:
        apply_move("unpoke", dd_unpoke_s1.value, dd_unpoke_s2.value)
    )
    btn_slide.on_click(lambda _:
        apply_move("slide", dd_slide_lms.value, tg_slide_om.value, tg_slide_ol.value, tg_slide_or.value)
    )
    btn_saddle.on_click(lambda _:
        apply_move("saddle", dd_saddle_s1.value, dd_saddle_s2.value)
    )
    btn_birth.on_click(lambda _:
        apply_move("birth")
    )
    btn_death.on_click(lambda _:
        apply_move("death", dd_death_loop.value) if dd_death_loop.value is not None else log("[death] no loop strands available")
    )

    def on_perm(_):
        try:
            d = ast.literal_eval(perm_text.value.strip())
            if not isinstance(d, dict):
                raise ValueError("Permutation must be a dict like {0:2, 2:0}")
            apply_move("permute", d)
        except Exception as e:
            log(f"[permute parse error] {e}")
    btn_perm.on_click(on_perm)

    frame.observe(lambda ch: render(), names="value")
    show_strand_graph.observe(lambda ch: render(), names="value")
    show_pd.observe(lambda ch: render(), names="value")
    plot_components_separately.observe(lambda ch: render(), names="value")

    refresh_dropdowns()
    render()

    header = W.VBox([
        W.HBox([play, frame, btn_jump_last]),
        W.HBox([show_strand_graph, show_pd, plot_components_separately]),
        W.HBox([pd_text, W.VBox([btn_load, btn_clear_log])]),
    ])

    ui = W.VBox([
        header,
        W.HBox([out_main, out_side]),
        tab,
        W.HTML("<b>Log</b>:"),
        out_log
    ])

    display(ui)
    return movie
