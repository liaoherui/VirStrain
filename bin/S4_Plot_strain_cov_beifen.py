import re
import os
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from _plotly_future_ import v4_subplots
import plotly.plotly as py
import  plotly.graph_objs as go
from plotly.subplots import make_subplots
#from plotly import tools

sfl=pd.read_csv('Mps_ps_depth_HIV.csv')
sf2=pd.read_csv('Ops_ps_depth_HIV.csv')
#print(sfl.columns)
# Get 
s=[]
for e in list(sfl.columns):
	if re.search('ID',e):continue
	name=re.split('_',e)[0]
	if name not in s:
		s.append(name)
so=[]
for e in list(sf2.columns):
	if re.search('ID',e):continue
	if re.search('None',e):continue
	name=re.split('_',e)[0]
	if name not in so:
		so.append(name)

#fig=plt.figure(figsize=(17,6))
#ax=fig.add_subplot(111)
#width=0.2

#sfl.plot('ID',s[0]+'_Freq',kind='bar')
#plt.xlabel("")
#plt.savefig('Test_cov.png')
#fig=py.iplot({'data':[Bar(x=sfl['ID'],y=sfl[s[0]+'_Freq'])],'layout':{'margin':{'b':300},'xaxis':{'tickangle':50}}})
#fig = tools.make_subplots(rows=2, cols=1, subplot_titles=['Depth', 'Label_Number'], shared_xaxes=True, shared_yaxes=False)
#traces=[go.Bar(x=sfl['ID'],y=sfl[s[0]+'_Freq']),go.Bar(x=sfl['ID'],y=sfl[s[0]+'_LNum'])]
# Extract 0 freq
i=0
zx=[]
zy=[]
for sw in sfl[s[0]+'_Freq']:
	if sw==0: 
		zx.append(sfl['ID'][i])
		zy.append(sw)
	i+=1
ft=open('Test.html','w+')
fig=go.Figure()
fig = make_subplots(specs=[[{"secondary_y": True}]])
fig.add_trace(go.Bar(x=sfl['ID'],y=sfl[s[0]+'_Freq'],name='depth'),secondary_y=False,)
#fig.add_trace(go.Scatter(x=zx,y=zy,mode='markers',opacity=0.5,marker=dict(size=5),name='Zero depth pos'),secondary_y=False,)
fig.add_trace(go.Scatter(x=sfl['ID'],y=(-1)*sfl[s[0]+'_LNum'],mode='lines+markers',opacity=0.5,marker=dict(size=5),name='(-1)*(strain number)'),secondary_y=True,)
fig.add_trace(go.Scatter(x=zx,y=zy,mode='markers',opacity=0.5,marker=dict(size=5,color='red'),name='Zero depth pos'),secondary_y=False,)
#fig2=go.Figure()
#fig2.add_trace(go.Bar(x=sfl['ID'],y=(sfl[s[0]+'_LNum'])))
'''
fig.append_trace(traces[0],1,1)
fig.append_trace(traces[1],2,1)
fig.layout.update(autosize=False,width=1000,height=500,)
mi=sfl[s[0]+'_LNum'].min()
ma=sfl[s[0]+'_LNum'].max()
print(mi)
fig['layout']['yaxis2'].update(title='',range=[mi,ma],autorange=False)
'''
ma=sfl[s[0]+'_LNum'].max()
ma2=sfl[s[0]+'_Freq'].max()
fig.layout.update(autosize=False,width=1200,height=300,title={'text':'Most Possible Strain: '+s[0],'xanchor': 'center'})
#fig['layout']['yaxis2']['autorange']
#fig['layout']['yaxis2']['autorange'] = "reversed"
fig.update_yaxes(range=[(-3)*ma,3*ma],secondary_y=True,)
#fig['layout']['yaxis2']['autorange'] = "reversed"
fig.update_yaxes(range=[(-1)*ma2,ma2],secondary_y=False,)
#with open('Test.html','a') as ft:
ft.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#ft.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#ft.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#ft.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#fig.write_html('Test.html')

