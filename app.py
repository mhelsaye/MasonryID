import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import io
import base64
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
from MasonryWalls_ID import draw_blocks, cross_section, solve_betaC,  calculate_point2, calculate_point3, calculate_pure_moment

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
# units 
mm = 0.001
m = 1
N = .001
kN = 1
kPa = 1

# Inputs 
H = 8 *m 
t = 240 *mm  
s = 1200 *mm 
db = 25 *mm
d = t/2 
fblock= 25 
faim = 0.6
fais = 0.85
emu = 0.003
k= 1

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Default values
H_default = 8000 *mm   # mm
t_options = np.array([140, 190, 240, 290]) *mm  # mm
t_default = 240 *mm   # mm
fblock_options = np.array([10, 15, 20, 25, 30])   # MPa
fblock_default = 25  # MPa
S_options = np.array([200,400,600,800,1000,1200, 1400, 1600]) *mm  # mm
S_default = 1200 *mm  # mm
bar_options = np.array([10, 15, 20, 25]) *mm  # mm
bar_default = 25 *mm  # mm
P_DL_default = 10 *kN # kN
P_LL_default = 0  *kN  # kN
e_default = 0  *mm # mm
W_default = 1 *kPa # kPa


# Layout
def generate_layout():

    return dbc.Container([

        html.H1("Masonry Wall Design Tool", className="text-center my-4", style={"font-size": "32px"}),

        dbc.Row([

            dbc.Col([

                html.H5("Inputs", className="text-primary", style={"font-size": "24px"}),

                dbc.Label("Wall Height[m]", style={"font-size": "18px"}),

                dbc.Input(id="input-H", type="number", value=H_default, step=1, min=0, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Wall Thickness [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-t", options=[{"label": f"{int(val*1000)} ", "value": val} for val in t_options], value=t_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Block Strength (fblock) [MPa]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-fblock", options=[{"label": f"{val} ", "value": val} for val in fblock_options], value=fblock_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Bar Spacing [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-S", options=[{"label": f"{int(val*1000)}", "value": val} for val in S_options], value=S_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Rebar Diameter [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-bar", options=[{"label": f"{int(val*1000)} ", "value": val} for val in bar_options], value=bar_default, style={"margin-bottom": "30px", "font-size": "16px"}),
                 html.Button("Check", id="btn-check", className="btn btn-primary mt-4", style={"font-size": "18px"}),

            ], width=4),

            dbc.Col([

                html.H5("Loads", className="text-primary", style={"font-size": "24px"}),

                dbc.Label("Dead Load  [kN]", style={"font-size": "18px"}),

                dbc.Input(id="input-P_DL", type="number", value=P_DL_default, step=0.1, min=0, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Live Load  [kN]", style={"font-size": "18px"}),

                dbc.Input(id="input-P_LL", type="number", value=P_LL_default, step=0.1, min=0, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Eccentricity  [mm]", style={"font-size": "18px"}),

                dbc.Input(id="input-e", type="number", value=e_default, step=0.1, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Wind Load [kPa]", style={"font-size": "18px"}),

                dbc.Input(id="input-W", type="number", value=W_default, step=0.1, style={"margin-bottom": "30px", "font-size": "18px"}),

               

            ], width=4),

            dbc.Col([

                html.Div(id="side-view", className="border", style={"height": "500px", "padding": "0", "overflow": "hidden", "backgroundColor": "transparent"}),

               html.Div(id="wall-image", className="border", style={"height": "100px", "padding": "0", "overflow": "hidden"}),

            ], width=4),

        ], className="my-4"),
                dbc.Row([

            dbc.Col([

                html.H5("Interaction Diagram", className="text-primary"),

                html.Div(id="interaction-diagram", className="border", style={"height": "600px"}),

            ], width=8),

        ]),


    ])

@app.callback(
    Output("wall-image", "children"),
    Input("dropdown-t", "value"),
    Input("dropdown-S", "value"),
    Input("dropdown-bar", "value")
)
def update_wall_image(t, s, bar):
    fig = draw_blocks(height=140, faceshell=30, width=400, num_blocks=5, gap=10, s=s/mm)
    
    ax = fig.gca()
    ax.set_aspect("auto")  # Automatically adjust aspect ratio
    ax.axis("off")  # Turn off axes completely
    
    # Adjust limits and padding to remove extra white space
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    buf = io.BytesIO()
    FigureCanvas(fig).print_png(buf)  # Save the figure as PNG
    data = base64.b64encode(buf.getvalue()).decode()
    
    return html.Img(src=f"data:image/png;base64,{data}", style={"width": "100%", "height": "100%", "object-fit": "contain"})


@app.callback(
    Output("side-view", "children"),
    [Input("input-H", "value"),
     Input("input-P_DL", "value"),
     Input("input-P_LL", "value"),
     Input("input-e", "value"),
     Input("input-W", "value")])

def update_side_view(H, P_DL, P_LL, e, W):
    # Default values if inputs are None
    if H is None:
        H = H_default
    if P_DL is None:
        P_DL = P_DL_default
    if P_LL is None:
        P_LL = P_LL_default
    if e is None:
        e = 0

    # Create the figure
    fig = Figure(figsize=(4, 8))
    ax = fig.add_subplot()
    ax.set_aspect("auto")  # Automatically adjust aspect ratio
    ax.set_xlim(-1, 2)  # Fixed x-axis range for consistent scaling
    ax.set_ylim(-0.5, 10.5)  # Fixed y-axis range for consistent scaling


    # Define the start points (x, y) of the arrows
    x_cor= -0.6
    x = np.array([x_cor, x_cor, x_cor, x_cor, x_cor,x_cor, x_cor, x_cor, x_cor])
    y = np.array([5,5.5,  6, 6.5 , 7 , 7.5, 8, 8.5, 9])

    # Define the direction and length of the arrows (u, v)
    Arr_L = 0.2
    Arr_D = 0
    u = np.array([Arr_L,Arr_L,Arr_L,Arr_L,Arr_L,Arr_L,Arr_L,Arr_L,Arr_L,Arr_L])
    v = np.array([Arr_D,Arr_D,Arr_D,Arr_D,Arr_D,Arr_D,Arr_D,Arr_D,Arr_D,Arr_D])

    # Add the arrows to the plot
    for i in range(len(x)):
        ax.arrow(x[i], y[i], u[i], v[i], head_width=0.1, head_length=0.1, fc="blue", ec="blue")

    # Coordinates of the triangle's vertices
    triangle_coords = [(-0.1,4.6), (0, 4.9), (0.1, 4.6)]

    # Create a Polygon patch
    triangle = patches.Polygon(triangle_coords, closed=True, color='black')
    ax.add_patch(triangle)
    # Draw the wall (fixed height)
    wall_thickness = 0.4  # Fixed wall thickness in plot units
    ax.add_patch(Rectangle((-wall_thickness / 2, 5), wall_thickness, 4, color="lightgrey"))

    # Dashed lines for height `h`
    ax.plot([0, 0], [5, 10], linestyle="dashed", color="black", linewidth=1.5)
    # ax.plot([-0.5, 1.5], [10, 10], linestyle="dashed", color="black", linewidth=1.5)
    # ax.plot([-0.5, 1.5], [5, 5], linestyle="dashed", color="black", linewidth=1.5)

    # Add height label
    ax.annotate(f"h = {H:.2f} m", xy=(0.4, 7), fontsize=12, color="black",
                rotation=90, va="center", ha="center", bbox=dict(facecolor="white", edgecolor="black", boxstyle="round"))
    
    # Add Wind label
    ax.annotate(f"W = {W:.2f} kPa", xy=(-0.7, 7), fontsize=12, color="black",
                rotation=90, va="center", ha="center")

    # Load arrows and labels
    arrow_length = 1  # Standard arrow length
    # Dead load arrow (P_DL)
    ax.arrow(e/1000, 10, 0, -arrow_length, head_width=0.2, head_length=0.3, fc="blue", ec="blue")
    ax.annotate(f"PDL = {P_DL:.2f} kN", xy=(e/1000+0.2, 9.8), fontsize=12, color="blue", va="center", ha="left")

    # Live load arrow (P_LL)
    # ax.arrow(0, 6, 0, arrow_length, head_width=0.2, head_length=0.3, fc="red", ec="red")
    ax.annotate(f"PLL = {P_LL:.2f} kN", xy=(e/1000+0.2, 9.3), fontsize=12, color="red", va="center", ha="left")

    # Eccentricity arrow (e)
    # ax.arrow(-0.3, 9.8, 0.6, 0, head_width=0.2, head_length=0.1, fc="purple", ec="purple")
    ax.annotate(f"e = {e:.2f} mm", xy=(0, 10.2), fontsize=12, color="purple", va="center", ha="center")

    # Turn off axis
    ax.axis("off")
    # Remove margins and set background to transparent
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)  # Remove margins
    ax = fig.gca()
    ax.axis('off')  # Remove axes
    ax.set_facecolor('none')  # Make axes background transparent
    fig.set_facecolor('none')  # Make figure background transparent

    # Save the figure to a buffer
    buf = io.BytesIO()
    FigureCanvas(fig).print_png(buf)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode()

    # Return as an image
    return html.Img(src=f"data:image/png;base64,{img_base64}", style={"width": "100%", "height": "auto"})


@app.callback(
    Output("interaction-diagram", "children"),
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-e", "value"),
    State("input-W", "value")
)



def update_interaction_diagram(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, e, W):
    if not n_clicks:
        return ""
        
    # Get cross section properties
    beff_m_1, beff_m_2, Aseff_m, bg_m, bug_m_1, bug_m_2, Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n ,E_m , ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf =cross_section(t, s,db)
    
    # Calculate maximum point
    PMax = solve_betaC(0.6,   fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t)
    betaC1 = float(PMax[0])
    
    # Calculate points for interaction diagram
    point2_results = calculate_point2(betaC1,faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=15)
    point3_results, Mr_y, Pr_y, ey = calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points=40)
    pure_moment = calculate_pure_moment(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t)
    
    
    # Create arrays for plotting
    M = [0] + [PMax[-1]] + [pt[4] for pt in point2_results] + [pt[4] for pt in point3_results] + [pure_moment[4]]
    P = [PMax[-2]] + [PMax[-2]] + [pt[3] for pt in point2_results] + [pt[3] for pt in point3_results] + [pure_moment[3]]

    lambda_h = k *H/t 
    # Loads Calculation 
    # Loads Calculation 
    e_P = np.maximum(0.1 * t, e)
    # Self Weight
    P_SW_mid = rho_SW  * H/2  # at mid span 
    P_SW_base = rho_SW  * H # at base

    # moment due to wind 
    M_lateral = W * H**2 / 8
    M_unf = (P_DL + P_LL) * e_P + M_lateral

    M_lateral_F1 = 0
    M_lateral_F2 = 0
    M_lateral_F3 = 0
    M_lateral_F4 = 1.4 *  M_lateral
    M_lateral_F5 = 1.4 *  M_lateral

    M_lateral_F = [M_lateral_F1, M_lateral_F2, M_lateral_F3, M_lateral_F4, M_lateral_F5]    

    # Factored Primary Moment at midspan 
    M_F1 = np.maximum(1.4 * P_DL * e_P /2, 1.4 * (P_DL  + P_SW_mid) * 0.1 *t)
    M_F2 = np.maximum(1.25 * P_DL * e_P /2+ 1.5 * P_LL * e_P /2, 1.25 * (P_DL  + P_SW_mid) * 0.1 *t + 1.5 * P_LL * 0.1 *t)
    M_F3 = np.maximum(0.9 * P_DL * e_P/2 + 1.5 * P_LL * e_P  /2, 0.9 * (P_DL  + P_SW_mid) * 0.1 *t + 1.5 * P_LL * 0.1 *t)
    M_F4 = np.maximum(1.25 * P_DL * e_P /2+ 1.4 * M_lateral, 1.25 * (P_DL  + P_SW_mid) * 0.1 *t)
    M_F5 = np.maximum(0.9 * P_DL * e_P /2+ 1.4 * M_lateral, 0.9 * (P_DL  + P_SW_mid) * 0.1 *t)

    M_F = [M_F1, M_F2, M_F3, M_F4, M_F5]

    # Factored Primary Axial Load @ midspan
    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid) 

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5]
    # top Moment 
    M1_F1 =  np.maximum(1.4 * P_DL * e_P, 1.4 * P_DL *0.1*t)
    M1_F2 =  np.maximum(1.25 * P_DL * e_P+ 1.5 * P_LL * e_P, 1.25 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
    M1_F3 =  np.maximum(0.9 * P_DL * e_P + 1.5 * P_LL * e_P, 0.9 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
    M1_F4 =  np.maximum(1.25 * P_DL * e_P , 1.25 * P_DL * 0.1 *t)
    M1_F5 =  np.maximum(0.9 * P_DL * e_P  , 0.9 * P_DL * 0.1 *t)

    Mtop_F= [M1_F1, M1_F2, M1_F3, M1_F4, M1_F5]
    # Bottom Moment 
    M2_F1 = np.maximum(0,  1.4 * (P_DL  + P_SW_base) *0.1*t)
    M2_F2 = np.maximum(0,  1.25 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
    M2_F3 = np.maximum(0,  0.9 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
    M2_F4 = np.maximum(0, 1.25 * (P_DL  + P_SW_base) * 0.1 *t)
    M2_F5 = np.maximum(0, 0.9 * (P_DL  + P_SW_base) * 0.1 *t)

    Mbot_F= [M2_F1, M2_F2, M2_F3, M2_F4, M2_F5]

    # Moment Ratio
    M_ratio =  np.minimum(np.array(Mtop_F), np.array(Mbot_F)) / np.maximum(np.array(Mtop_F), np.array(Mbot_F))

    ev = np.array([M_F1, M_F2, M_F3, M_F4, M_F5]) / np.array([P_F1, P_F2, P_F3, P_F4, P_F5])
    # Rigidty Coefficient 
    betad = P_DL * e_P / M_unf
    betad=0
    faiE = 0.75
    EI_eff_raw = E_m * (0.25 * I_gross_eff - (0.25 * I_gross_eff - I_cr_eff)* ((ev - ek) / (2 * ek))) 
    EI_eff = np.clip(EI_eff_raw, E_m * I_cr_eff, 0.25 * E_m * I_gross_eff)
    Rigidty_c = faiE * EI_eff / (1+0.5*betad)

    Pcr = np.pi**2 * Rigidty_c/ (k*H)**2*1000
    # Lateral Force Coefficient
    if lambda_h < 30:
        Cm = np.where(np.array(M_lateral_F) / np.array(M_F) > 0.5,
                    1,
                    np.maximum(0.6 + 0.4 * M_ratio, 0.4))
    else: 
        Cm = [1] * len(M_F)

    MagFactor=  np.array(Cm) / (1-(np.array(P_F)/np.array(Pcr)))

    # Design Moment

    Mt_F = np.array(M_F) * np.array(MagFactor)

    # 
    Mt_F_list = [float(val) for val in Mt_F]
    M_F_list = [float(val) for val in M_F]
    P_F_list = [float(val) for val in P_F]
    
    # Create the figure
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=M, y=P, mode="lines", name="Interaction Curve"))

    fig.add_trace(go.Scatter(x=Mt_F_list, y=P_F_list, mode="markers", name="Total Moments"))
    fig.add_trace(go.Scatter(x=M_F_list, y=P_F_list, mode="markers",marker=dict(color='blue',symbol='x'), name="Primary Moment"))
    fig.update_layout(
        xaxis_title="Moment (kN⋅m)",
        yaxis_title="Axial Force (kN)",plot_bgcolor="white",
        xaxis=dict(gridcolor="lightgrey", range=[0, 1.1*max(max(M),max(Mt_F_list))],    
                   zeroline=True, zerolinecolor="black",zerolinewidth=1),
        yaxis=dict(gridcolor="lightgrey", range=[0, 1.1*max(max(P),max(P_F_list))],
                   zeroline=True, zerolinecolor="black",zerolinewidth=1),
        showlegend=True,
        margin=dict(l=0, r=0, b=0, t=0.1),  # Remove margins
                      paper_bgcolor='rgba(0,0,0,0)'
    )
    
    return dcc.Graph(figure=fig)

app.layout = generate_layout

if __name__ == "__main__":
    app.run_server(debug=True)
