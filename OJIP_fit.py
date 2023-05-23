import plotly.graph_objects as go
import pandas as pd
import numpy as np
import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import base64
import io
from scipy import optimize
import dash_bootstrap_components as dbc

from dash import dash_table
from mvgavg import mvgavg



x = np.linspace(0, 10, 50)
y = np.exp(-x)

df_init = pd.DataFrame(np.array([x,y]).T, columns = ["time", "fluorescence"])
df_init=df_init.to_dict()



def pre_process(time_array, fluo, N_mvg = 10, N_log = 1000 ):
    #remove blank
    blank = np.mean(fluo[0:10])
    fluo = fluo-blank

    #stop at max
    #stop_fluo = int(np.argmax(fluo))
    #fluo = fluo[0:stop_fluo]
    #time_array = time_array[0:stop_fluo]

    #normalise
    fluo_ref = fluo/fluo.max()    

    #collect the fluorescence rise: during the jump of light
    ind_ref = (fluo_ref>0.1)
    #perform moving average
    binit = True
    t = mvgavg(time_array[ind_ref], N_mvg, binning = binit)
    y = mvgavg(fluo[ind_ref], N_mvg, binning = binit)
    
    # start at 0
    t-=t[0]
    
    #logarithmic subsampling to accelerate the fit
    ind= np.unique(np.logspace(np.log10(1), np.log10(len(t)-1), N_log).astype(int))
    t = t[ind]
    y = y[ind]

    return t, y


def sigmoidal_OJIP(parameters, tdata):
    F0 = parameters[0]
    Aoj = parameters[1]
    koj = parameters[2]
    soj = parameters[3]
    Aji = parameters[4]
    kji = parameters[5]
    sji = parameters[6]
    Aip = parameters[7]
    kip = parameters[8]
    sip = parameters[9]
    y = F0 + Aoj*(1-np.exp(-koj*tdata))**soj + Aji*(1-np.exp(-kji*tdata))**sji +  Aip*(1-np.exp(-kip*tdata))**sip
    
    return y
     
def multiexp_fit(t, y):
    """ triexponential sigmoidal fit of the fluorescence rise based on Joly & Carpentier 2009"""

    #initial parameters based on Joly & Carpentier, 2009
    dF = y.max()-y.min()
    x0 = [y.min(), dF/2, 5E3, 1.24, dF/4, 0.06E3, 1.2, dF/4, 0.0023E3, 8.2]

    parameters_estimated = optimize.least_squares(residuals,  x0, bounds = ([-1e5,-1e5, 0, 1,-1e5, 0, 1,-1e5, 0, 1], [1e5,1e5, 1e5, 20,1e5, 1e5, 20,1e5, 1e5, 20]),
                                args = (t, y, sigmoidal_OJIP))
    
    #recover the characteristic time of the first phase (O-J)
    tau = 1/parameters_estimated.x[2]
    
    ypred = sigmoidal_OJIP(parameters_estimated.x, t)
    
    return tau, ypred
    


def get_fit(t, y):    
    time_spread = t.max()-t.min()
    start = np.mean(y[0])
    stop = np.mean(y[-10:])
    x0 = [start, 1/time_spread, stop]
    t = t-t[0]

    parameters_estimated = optimize.least_squares(residuals,  x0, bounds = (-1e8,1e8),
                                args = (t, y, exp_decay))

    tau = parameters_estimated.x[1]

    pos_tau = find_nearest(t, tau)

    x0 = parameters_estimated.x #initial guess: parameters from previous fit
    #second fit
    parameters_estimated  = optimize.least_squares(residuals,  x0, bounds = (-1e9,1e9),
                                            args = (t[0:int(pos_tau*5)], y[0: int(pos_tau*5)], exp_decay))

    return  parameters_estimated.x


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



    
def exp_decay(parameters, xdata):
    '''
    Calculate an exponential decay of the form:
    S = a * exp(-xdata/b)
    '''
    A = parameters[0]
    tau = parameters[1]
    y0 = parameters[2]
    return A * (1 - np.exp(-xdata/tau))**1.24 + y0


def residuals(parameters, x_data, y_observed, func):
    '''
    Compute residuals of y_predicted - y_observed
    where:
    y_predicted = func(parameters,x_data)
    '''
    return func(parameters,x_data) - y_observed



# Define a function to calculate the value based on the selected component
def calculate_value(sigma, params):
    # Replace this with your own calculation logic based on the chemical and wavelength
    return 1e6/(sigma*params[1])

# Create the Dash app instance
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define the layout of the app

controls = dbc.Card(
    [
        html.Div(
            [

        html.Div(id='select-wavelength', 
                                children=[html.Div('Select the excitation wavelength used:', style={'font-weight': 'bold', 'margin-right': '24px',
                                'font-size': '24px'})]),
        dcc.Input(
            id='wavelength-dropdown',
            type='number',
            inputMode='numeric',
            placeholder='Enter an integer',
            value = 470,
                ),

        html.Div(id='output-container', 
                                children=[html.Div('Sigma (m²/mol):', style={'font-weight': 'bold', 'margin-right': '10px', "color":"darkred"})]),
        html.Div(id='sigma-value', 
                 style={'display': 'inline-block', 'font-size': '24px', 'vertical-align': 'middle', "color":"darkred"},
                 ),
            ]),
        
         html.Div(id='direct-table', 
                                children=[html.Div(html.Div([
                                            html.P([
                                                html.Strong('Upload your table:')
                                            ], style={'font-size':'24px'}),
                                            html.Ul([
                                                html.Li('.csv file'),
                                                html.Li('column separator: comma ","'),
                                                html.Li('decimal separator: dot "."')],
                                             style={'display': 'inline-block', 'font-size': '18px', 'vertical-align': 'middle'}),
                                            ]),
                                            
                                            
                                
                                )]), 
        

        dcc.Upload(         id='upload-data',
                            children=html.Div([
                                     'Drag and Drop or ',
                            html.A('Select Files')
                        ]),
                            style={
                                'width': '80%',
                                'height': '60px',
                                'lineHeight': '60px',
                                'borderWidth': '1px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '10px',
                            },
                            multiple=False,    
                        ),
           

    html.Div(id='output-container3', 
                    children=[html.Div('Time constant (s):', style={'font-weight': 'bold', 'margin-right': '10px'})]),

    html.Div(id='tau-value', 
                
                style={'display': 'inline-block', 'font-size': '24px', 'vertical-align': 'middle'}),

        html.Div(id='output-container2', 
                    children=[html.Div('Light intensity (µE/m²/s):', style={'font-weight': 'bold', 'margin-right': '10px'})]),

    html.Div(id='intensity-value-eins', 
                style={'display': 'inline-block', 'font-size': '24px', 'vertical-align': 'middle'},
                ),
    html.Div(id='output-container5', 
                    children=[html.Div('Light intensity (mW/mm²):', style={'font-weight': 'bold', 'margin-right': '10px'})]),

    html.Div(id='intensity-value-watt', 
                style={'display': 'inline-block', 'font-size': '24px', 'vertical-align': 'middle'},
                )  
    ]   
)
            

loading =  dbc.Card(
    [   
        dcc.Loading(
                id="loading",
                type="cube",
                children=[dcc.Store(id='data-store', data = df_init),           
            ],
                style={
                                'width': '80%',
                                'height': '60px',
                            }
        ),

        dcc.Loading(
                id="loading2",
                type="graph",
                children=[dcc.Store(id='fit-store', data = "None")],
                style={
                                'width': '80%',
                                'height': '60px',
                            },
            
        ),
    ]
    )


upload_table = dbc.Card(
    [   
        html.Div(id='output-table')

            ]
        )





graph =  html.Div([

            html.Div(id='select-axis', 
                                children=[html.Div('Select the X and Y column names:', style={'font-weight': 'bold', 'margin-right': '10px',
                                "font-size":"24px"})]),

            dcc.Dropdown(
                id='x-axis-dropdown',
                options=[],
                value=None,
                placeholder="Select X-axis Column",
                style={'margin': '10px'}
                    ),
            dcc.Dropdown(
                id='y-axis-dropdown',
                options=[],
                value=None,
                placeholder="Select Y-axis Column",
                style={'margin': '10px'}
            ),

            html.Div(id='select-smooth', 
                                children=[html.Div('Select the smoothing window size and logarithmic subsampling:', style={'font-weight': 'bold', 'margin-right': '10px',
                                "font-size":"24px"})]),
            dcc.Input(
                id='smooth-dropdown',
                type='number',
                inputMode='numeric',
                placeholder='Enter an integer',
                value = 10,
                ),

            dcc.Input(
                id='log-dropdown',
                type='number',
                inputMode='numeric',
                placeholder='Enter an integer',
                value = 10000,
                ),

    dcc.Graph(id='data-plot',
                    style={
                    #'height': '500px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                    },
                
                    ),

            ])   



app.layout = dbc.Container(
    [
        html.H1("Light calibration"),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(controls, md=4),
                dbc.Col(graph, md=4),
                dbc.Col(loading, md=4)
            ]),

               
      ],
            fluid=True,

)




"""SCALAR VALUES"""
# Define the callback function that updates the output value based on the selected chemical and wavelength
@app.callback(
    Output('sigma-value', 'children'),
    Input('wavelength-dropdown', 'value')
)
def update_output_value(wl):
    if wl is None:
        return ""
    else:
        return '{:.1e}'.format(sigma_spectra[wavelength==wl][0])
    
# Define the callback function that updates the output value based on the selected chemical and wavelength
@app.callback(
    Output('tau-value', 'children'),
    Input('fit-store', 'data')
)
def update_tau_value(params):
    if params is None:
        return ""
    else:
        return '{:.1e}'.format(params['params_exp'][1])
        
@app.callback(
    Output('intensity-value-watt', 'children'),
    Input('fit-store', 'data'),
    Input('wavelength-dropdown', 'value'),

)
def update_output_value(dico, wl):
    if wl is None  or dico is None:
        return ""
    else:
        sigma = sigma_spectra[wavelength==wl]
        value = calculate_value( sigma, dico['params_exp'])[0]
        return '{:.1e}'.format(value/1000*120/int(wl))
    
@app.callback(
    Output('intensity-value-eins', 'children'),
    Input('fit-store', 'data'),
    Input('wavelength-dropdown', 'value')
)
def update_output_value(dico, wl):
    if wl is None or dico is None:
        return ""
    else:
        sigma = sigma_spectra[wavelength==wl]
        value = calculate_value(sigma, dico['params_exp'])[0]
        return '{:.1e}'.format(value)

"""DATA STORAGE"""
@app.callback(
    Output('data-store', 'data'),
    Input('upload-data', 'contents'),
    Input('upload-data', 'filename')
)
def update_storage(contents, filename):
    if contents is not None:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            if 'csv' in filename:
                # Assume that the user uploaded a CSV file
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))

            elif 'xls' in filename:
                # Assume that the user uploaded an Excel file
                df = pd.read_excel(decoded, engine = 'openpyxl', encoding='ISO-8859-1')

            elif 'txt' or 'tsv' in filename:
                # Assume that the user upl, delimiter = '\t'
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), delimiter = '\t')
            
            else:
                return html.Div([
                    'The uploaded file format is not supported. Please upload a CSV, Excel, TXT or TSV file.'
                ])
            


            return df.to_dict()
        
        except Exception as e:
            print(e)
            return html.Div([
                'There was an error processing this file.'
            ])


def read_table(df):
        df = pd.DataFrame(df)
        return html.Div([
            dash_table.DataTable(
                data=df.to_dict('records'),

                columns=[{'name': i, 'id': i} for i in df.columns],
                        
                style_table={'overflowX': 'scroll', 
                             'maxHeight': '200px',
                             'maxWidth': '100%',
                             'overflowY': 'scroll'},

            )])
"""
@app.callback(
    Output('output-table', 'children'),
    Input('data-store', 'data'),
)
def update_table(df):
        if df is not None:
            children = [read_table(df)]
            return children
"""
            
@app.callback(
    dash.dependencies.Output('x-axis-dropdown', 'options'),
    dash.dependencies.Output('y-axis-dropdown', 'options'),
    dash.dependencies.Input('data-store', 'data')
)
def update_dropdowns(df):
    if df is not None:
        df = pd.DataFrame(df)
        x_axis_options = [{'label': col, 'value': col} for col in df.columns]
        y_axis_options = [{'label': col, 'value': col} for col in df.columns]
        return x_axis_options, y_axis_options
    else:
        return [], []

@app.callback(
    Output('fit-store', 'data'),
    Input('data-store',"data"),
    Input('x-axis-dropdown', 'value'),
    Input('y-axis-dropdown', 'value'),
    Input('smooth-dropdown', 'value'),
    Input('log-dropdown', 'value'),

)
def update_fit(df, key_time, key_fluo, N_mvg, N_log):
    if key_time is None or df is None:
        return None
    else:
        df = pd.DataFrame(df)

        t, y = pre_process(df[key_time], df[key_fluo], N_mvg = N_mvg, N_log = N_log)
        tau, ypred =  multiexp_fit(t,y)
        print(tau)
        pos = 3*tau
        pos_tau = find_nearest(t, pos)
        params = get_fit(t[:pos_tau], y[:pos_tau])
        print(params[1])
        return {"y_JC":ypred, "params_exp":params, "t": t, "y": y}


# Define the callback function that collects the file and reads it
@app.callback(
    Output('data-plot', 'figure'),
    Input('fit-store',"data"),
    Input('x-axis-dropdown', 'value'),
    Input('y-axis-dropdown', 'value'),

)
def update_figure(dico, key_time, key_fluo):
            fig = go.Figure()

            fig.update_layout(
                    plot_bgcolor='white',  # Set the plot background color to white

                    xaxis=dict(showgrid=False),  # Hide the x-axis grid lines
                    yaxis=dict(showgrid=False),  # Hide the y-axis grid lines
)
            if dico is None or key_time is None or key_fluo is None:
                return fig
            else:
                X = dico['t']
                Y = dico['y']

                fig.add_trace(go.Scatter(
                    
                    x=X, 
                    y=Y, 
                    name = "raw",
                    mode = "markers",
                    marker_color='rgba(152, 0, 0, .8)',
                    )
                )
                x_pred = np.logspace(np.log10(np.min(X)), np.log10(np.max(X)), 1024)
                y_pred = exp_decay(dico["params_exp"], x_pred)
                fig.add_trace(go.Scatter(
                                x = x_pred.tolist(), 
                                y=y_pred.tolist(), 
                                name = "exponential fit",
                                mode = "lines",
                                line_color = "rgba(0,0,0,0.6)",
                                )
                )
                fig.add_trace(go.Scatter(
                                x = X, 
                                y=dico['y_JC'], 
                                name = "JC fit",
                                mode = "lines",
                                line_color = "rgba(1,0,0,1)",
                                )
                )
                fig.update_layout(
                    xaxis_title=key_time,  # Set the x-axis label
                    yaxis_title=key_fluo,  # Set the y-axis label
                    xaxis_type='log'
                )
                
                return fig


sigma_spectra = np.array([1080320.36324794, 1096525.82659956, 1112731.28995119,
       1128936.75330282, 1145142.21665445, 1161347.68000608,
       1177553.1433577 , 1193758.60670933, 1209964.07006095,
       1226169.53341258, 1242374.99676421, 1258580.46011584,
       1274785.92346746, 1290991.38681909, 1307196.85017072,
       1323402.31352234, 1339607.77687397, 1373153.35502663,
       1411905.79143144, 1451673.70138916, 1478670.56466059,
       1515274.3414749 , 1548152.43008973, 1539517.61161135,
       1583957.42642031, 1613578.00644628, 1643833.38937371,
       1665504.34583023, 1691833.2616521 , 1713721.44702742,
       1734676.01453642, 1757146.36611417, 1765487.82513961,
       1767094.41048528, 1783768.76508012, 1789866.41365457,
       1791637.52973493, 1801831.56652883, 1815463.26548964,
       1845030.58414562, 1846578.66360136, 1863979.89956837,
       1869297.85274318, 1902973.70718976, 1922833.10302972,
       1932393.08574838, 1920830.76990208, 1929339.3557416 ,
       1943226.70453895, 1946433.07570567, 1948326.3226245 ,
       1943552.35931765, 1929795.11622905, 1934624.02574136,
       1916399.22013779, 1882641.23395522, 1857959.33610857,
       1806948.59035621, 1753736.01031018, 1718573.34900475,
       1680012.77893194, 1631544.86418464, 1590600.18960301,
       1553916.54869223, 1517889.90183009, 1482708.31324421,
       1450548.18458663, 1429028.00325783, 1397699.82415812,
       1374989.16295397, 1374790.54338109, 1367168.5460158 ,
       1368033.77665265, 1379182.02567163, 1405920.42837571,
       1413738.91475073, 1434123.67956377, 1448988.46172296,
       1473030.9927783 , 1492334.03517098, 1497323.22892746,
       1503626.64754565, 1515327.58272851, 1523428.05879884,
       1529298.19431568, 1542857.14285714, 1542576.79188756,
       1541548.89922265, 1550567.31610827, 1536215.83950724,
       1512423.57573875, 1499206.69943644, 1485409.30531317,
       1479630.06861839, 1470755.5418759 , 1462681.65170448,
       1458408.4264974 , 1448040.64278621, 1438886.96113299,
       1420759.70904844, 1398294.20753932, 1372085.2049734 ,
       1333710.62395713, 1314941.31821197, 1277702.45173747,
       1257359.1401519 , 1220428.06595823, 1186748.58146191,
       1132252.80887725, 1091719.85740871, 1044718.58630797,
        996027.48626862,  951288.91747806,  904267.35425221,
        867574.78031876,  819640.33542742,  782829.37362567,
        745628.94470458,  715579.50112647,  686670.54770786,
        656777.46277505,  630307.47555064,  598427.46076374,
        566948.8898423 ,  543721.49269629,  521713.03584591,
        498082.13220402,  478332.80220379,  465021.58937712,
        445962.2678377 ,  425367.69197932,  409067.69657186,
        388986.76297839,  371871.3247358 ,  362868.2130055 ,
        352515.7156202 ,  340674.78530855,  334304.10521768,
        324257.76779525,  318966.60588279,  312956.6737609 ,
        309455.17324309,  309357.25769934,  305914.35345918,
        303200.38196645,  298094.6168562 ,  297644.82081274,
        292007.94069656,  297522.06926282,  295377.64262062,
        297085.45001645,  290429.1900518 ,  290534.31619159,
        291988.25802339,  291361.52979218,  294164.43232791,
        297278.81265898,  296679.11267928,  298932.41951389,
        301271.97836244,  302210.71624825,  305234.91677873,
        303828.20393297,  303788.05741674,  302627.86066651,
        304607.7980444 ,  299006.6492145 ,  302964.0241715 ,
        295382.91699819,  294406.87787231,  297872.26589567,
        300546.77447763,  303939.46652755,  303214.59489387,
        307395.09672919,  309366.33209283,  308924.85321731,
        313385.53114165,  317841.0477302 ,  325406.87857612,
        335877.30969035,  337368.2387927 ,  355642.58734382,
        358300.44169331,  367389.7945327 ,  375952.1663928 ,
        371255.72362331,  373648.53436215,  381892.83691055,
        380134.28410989,  387538.58267784,  392092.97517867,
        397844.23723147,  401655.75777558,  402206.89873785,
        405823.47024202,  404579.67076031,  402660.25674658,
        400331.65750322,  409568.91526944,  409465.29261404,
        417096.47974219,  417025.67367902,  421364.50223237,
        432240.59058115,  440388.65891419,  438779.67464752,
        440600.06346725,  446711.77753169,  456313.13167753,
        463354.87973147,  471658.83493694,  471477.21567323,
        482308.17816221,  486566.28537259,  489923.87621808,
        491983.72936474,  500368.59753164,  502571.00243522,
        505745.62455663,  503845.02961724,  506239.80310324,
        503651.24769957,  515835.76797197,  508511.24269749,
        515908.66986217,  508934.37201447,  513040.91397102,
        515390.89931683,  516475.20390449,  510321.36941409,
        512654.37736637,  509188.73522681,  508916.862424  ,
        514044.00659683,  513211.21717674,  520218.02740401,
        523064.84317152,  525364.24563228,  520416.34924037,
        525419.9910144 ,  527250.57854259,  527839.2979549 ,
        524666.47266005,  524180.28722176,  531624.19461285,
        527817.78088701,  534094.5718205 ,  539457.40050247,
        545880.14718371,  552107.79608908,  550494.33419321,
        557452.63669133,  580712.20522034,  605173.29073261,
        625164.06975426,  648216.08108103,  676233.06808596,
        691287.70760725,  735439.29711113,  766771.55370299,
        812219.08118596,  849139.64705903,  865441.92866851,
        892126.74625705,  917447.6092084 ,  936198.20234079,
        958605.81126973,  974184.78761212,  989773.78017848,
       1001107.03224672, 1016398.75442679, 1031686.71451576,
       1049660.39254716, 1073591.0810498 , 1055642.37679175,
       1134054.40309067, 1158380.66171228, 1152754.41263563,
       1189887.1017159 , 1211132.09693735, 1240335.41347613,
       1276705.44341248, 1283540.40693822, 1309489.70098349,
       1323724.42575004, 1326903.60188728, 1339888.60744422,
       1320294.97796485, 1329808.07135887, 1316872.04570731])

wavelength = np.array([385., 386., 387., 388., 389., 390., 391., 392., 393., 394., 395.,
       396., 397., 398., 399., 400., 401., 402., 403., 404., 405., 406.,
       407., 408., 409., 410., 411., 412., 413., 414., 415., 416., 417.,
       418., 419., 420., 421., 422., 423., 424., 425., 426., 427., 428.,
       429., 430., 431., 432., 433., 434., 435., 436., 437., 438., 439.,
       440., 441., 442., 443., 444., 445., 446., 447., 448., 449., 450.,
       451., 452., 453., 454., 455., 456., 457., 458., 459., 460., 461.,
       462., 463., 464., 465., 466., 467., 468., 469., 470., 471., 472.,
       473., 474., 475., 476., 477., 478., 479., 480., 481., 482., 483.,
       484., 485., 486., 487., 488., 489., 490., 491., 492., 493., 494.,
       495., 496., 497., 498., 499., 500., 501., 502., 503., 504., 505.,
       506., 507., 508., 509., 510., 511., 512., 513., 514., 515., 516.,
       517., 518., 519., 520., 521., 522., 523., 524., 525., 526., 527.,
       528., 529., 530., 531., 532., 533., 534., 535., 536., 537., 538.,
       539., 540., 541., 542., 543., 544., 545., 546., 547., 548., 549.,
       550., 551., 552., 553., 554., 555., 556., 557., 558., 559., 560.,
       561., 562., 563., 564., 565., 566., 567., 568., 569., 570., 571.,
       572., 573., 574., 575., 576., 577., 578., 579., 580., 581., 582.,
       583., 584., 585., 586., 587., 588., 589., 590., 591., 592., 593.,
       594., 595., 596., 597., 598., 599., 600., 601., 602., 603., 604.,
       605., 606., 607., 608., 609., 610., 611., 612., 613., 614., 615.,
       616., 617., 618., 619., 620., 621., 622., 623., 624., 625., 626.,
       627., 628., 629., 630., 631., 632., 633., 634., 635., 636., 637.,
       638., 639., 640., 641., 642., 643., 644., 645., 646., 647., 648.,
       649., 650., 651., 652., 653., 654., 655., 656., 657., 658., 659.,
       660., 661., 662., 663., 664., 665., 666., 667., 668., 669., 670.,
       671., 672., 673., 674., 675.])

# Run the app
if __name__ == '__main__':
    app.run_server(debug=False)
