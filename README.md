# A2TEA WebApp

<img src="a2tea_hexsticker.png" align="right" alt="" width="120" />

**Current setup instructions:**
  
A working conda and preferentially mamba installation (in conda base environment) is prerequisite:
  
Create a new conda environment & once activate it:  
```
conda create --name phaunos  
conda activate phaunos
```

Install necessary packages - R version 4.1.3!
```
mamba install -c conda-forge -c r r-base==4.1.3 jupyterlab r-irkernel
#also these - these required packages for ggmsa currently fail during R installation...
mamba install -c conda-forge r-gert r-proj4 r-devtools
```

Clone and navigate to the A2TEA WebApp directory
  

### Working on a local machine:
```
jupyter lab
```


### If you want to work on a remote machine - small guide to port forwarding (cs03 = alias for remote server):  
For this to work you need:

- install jupyter(lab) on both remote and local machine - e.g. `mamba install -c anaconda jupyter`
- on local machine install net-tools `sudo apt-get install net-tools`
- to make sure also install programming language on both local & remote machine (incl. the language specific jupyter kernel) - `conda install -c r r-irkernel`


On **remote machine** start jupyter and add the following options.  
If port 8888 isn't usable another one will chosen automatically (take notice of the number! - e.g. 8889).
```
jupyter lab --no-browser --port=8888
```

You will be greeted by a similar message to this:
```
 To access the server, open this file in a browser:
        file:///data/home/stoecker/.local/share/jupyter/runtime/jpserver-2687914-open.html
    Or copy and paste one of these URLs:
        http://localhost:8888/lab?token=7a4bf8eff1b9ca48376d61e6a2bc23179c39f565c21276db
     or http://127.0.0.1:8888/lab?token=7a4bf8eff1b9ca48376d61e6a2bc23179c39f565c21276db
```

Enter the following command on your **local machine**.
**Adapt port 8888 to the port displayed in the message on the remote machine!**

```
ssh -N -f -L localhost:8888:localhost:8888 cs03

```

If you stuck to the steps you can simply click (e.g. right-click -> open in browser) on one of the last two links given by the remote machine (these include the security token) and you will open jupyter lab on your local machine ;D  


Open the shiny_tests.ipnyb and install all required R packages (if first time); load packages; then run both ui and server functions.  
The last important cell actually runs the App.  
```
#Set port to the shiny specific port - NOT the port jupyter uses!
phaunos <- shinyApp(ui, server)
runApp(phaunos, 
       #host = '131.220.109.99', 
       port=8891)
```
We need a **second** port forward - so in this case we would need to type in our local machine `ssh -N -f -L localhost:8891:localhost:8891 cs03`.  
  
If everything worked you will be greeted with no errors but a blue link that you can click and which will open the WebApp on your local machine ;D
  
  
### Finishing work and closing connections
Remote/Server: It's enough to simply close the jupyter notebook and WebApp  
Local machine: 
Check open connections per port (you might need to use sudo)  
`netstat -lpn |grep :8891` will show you the Process ID of the tunnel connection, e.g.: 406659/ssh   
Kill connection with `kill 406659`  


### Further notes:  
- If sth. doesn't seem to work with the port, remember that perhaps are people are using them or you yourself might still have connections running!

