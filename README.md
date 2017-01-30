<ul>
 <li> ‘coreHMR.xml’ is the SBML for the metabolic model used in Damiani et al. </li>
 <li> ‘popFBA.m’ is the MATLAB function to perform Linear Optimisation on a population model composed of Npops clones of a given SBML model as described in Damiani et al. </li>
 <li> ‘createPopModel.m’ is the function called by popFBA that returns a population model structure from a metabolic model structure</li>
</ul>

<strong>coreHMR.xml</strong><br />
————————————————————————————————————————————————<br />
<p>
The exchange reactions are in the form of: metabolite[s] &lt;=&gt; ''<br>
The uptake/secretion reactions are in the form of: metabolite_in[c] &lt;=&gt; 
metabolite_out[s]<br />
The internal reactions are in the form of: reagent[compartment] &lt;=&gt; 
product[compartment]<br />
<br />
<br />
<strong>popFBA</strong><br />
</p>
————————————————————————————————————————————————<br />
[popModel, singleModel, optFlux] = popFBA(nameSBML, nameExRxns, nameCoopRxn, 
nPop, CharExtComp, otherFeat, rxnsFeat, metsFeat)<br />
<br />
popFBA returns a Matlab structure for the input SBML model (singleModel) and for 
the created population model (popModel) and the solution of the FBA problem for 
the popModel.<br />
see Matlab help for more details on function arguments.<br />
<br />
<br />
<em>USAGE EXAMPLE</em>:<br />
[modelPop, modelSingle, fluxese] = popFBA('namefile.xml', nameExRxns, 
nameCoopRxn, 's', 10, &lt;otherFeat&gt;, &lt;rxnsFeat&gt;, &lt;metsFeat&gt;)<br />
<br />
where<br />
<br />
nameExRxns = { ...<br />
'Ex_biomass',...<br />
'Ex_glucose',...<br />
'Ex_O2',...<br />
'Ex_glutamine',...<br />
'Ex_folate',...<br />
'Ex_lactate',...<br />
'Ex_urea',...<br />
'Ex_glutamate',...<br />
'Ex_NH3',...<br />
'Ex_H2O',...<br />
'Ex_Hs',...<br />
'Ex_putrescines',...<br />
'Ex_arginines',...<br />
'Ex_guanidinoacetatec',...<br />
'Ex_mercaptopyruvatec',...<br />
'Ex_phenylpyruvatec',...<br />
'Ex_4methylthio2oxobutanoicacidc',...<br />
'Ex_4methyl2oxopentanoatec',...<br />
'Ex_2oxo3methylvaleratec',...<br />
'Ex_citratec',...<br />
'Ex_CO2c',...<br />
'Ex_Pi',...<br />
'Ex_10formylTHFm'};<br />
<br />
nameCoopRxn = { ...<br />
'UptakeGlucose',...<br />
'LactateExchange',...<br />
'Ex_H',...<br />
'UptakeOxygen',...<br />
'AQP6',...<br />
'TransportCyNH3',...<br />
'UptakeGlutamine',...<br />
'DemandGlutamate',...<br />
'UptakeGlutamate',...<br />
'UptakeArginine',...<br />
'Ex_putrescine',...<br />
'SLC14A2',...<br />
'UptakeFolate',...<br />
'biomass_synthesis',...<br />
'GCK_HKDC1_HK1_ADPGK_HK2_HK3'};<br />
<br />
and 's' identify the tumor environmental compartment<br />
<br />
The set nameExRxns will be unchanged.<br />
<br />
The set nameCoopRxn will be cloned n = Npops times, as in the examples below:
<br />
metabolite_in_0[c] &lt;=&gt; metabolite_out[s], <br />
metabolite_in_1[c] &lt;=&gt; metabolite_out[s], <br />
..., <br />
metabolite_in_(n-1)[c] &lt;=&gt; metabolite_out[s]<br />
<br />
all remaining internal reactions will be cloned n = Npops time, as in the 
examples below: <br />
reagent_0[compartment] &lt;=&gt; product_0[compartment],<br />
reagent_1[compartment] &lt;=&gt; product_1[compartment],<br />
...,<br />
reagent_(n-1)[compartment] &lt;=&gt; product_(n-1)[compartment].<br />
<br />
If an external micro environment (s) is not already included in the single 
model, the cooperation reactions only should be indicated (no Exchange 
reactions).<br />
The algorithm will create exchange reactions for the metabolites involved in the 
cooperation reactions as in the example below: <br />
from metabolite_in[c] &lt;=&gt; '' to metabolite_in_0[c] &lt;=&gt; 
metabolite_out[s];metabolite_in_1[c] &lt;=&gt; metabolite_out[s] and so on; and 
metabolite[s] &lt;=&gt; ''<br />
<br />
<em>OPTIONAL ARGUMENTS</em><br />
otherFeat = {'genes','modelVersion'};<br />
rxnsFeat = {'rev','c','lb','ub','rules','grRules','rxnNames'};<br />
metsFeat = {'metNames','metCharge','metCHEBIID','metKEGGID','b'};<br />
<br />
Only the fields indicated will be cloned in the popModel<br />
