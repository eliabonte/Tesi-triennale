# "Analisi di prestazioni di comfort di diverse strategie di controllo di sospensione semi-attiva"

# Sommario
Le sospensioni sono una delle parti fondamentali di un veicolo, costituiscono il collegamento principale tra la strada e il telaio e sono responsabili della sicurezza e del comfort del passeggero a bordo.
In questo lavoro sono state analizzate, tramite simulazioni su un modello quarter car, le prestazioni di comfort di diverse strategie di controllo per le sospensioni semi-attive di tipo magnetoreologico.
La taratura dei parametri delle leggi di controllo è stata fatta tramite l’algoritmo di ottimizzazione bayesiana.
L’obiettivo iniziale era quello di osservare le prestazioni di due diverse leggi di controllo, lineare e quadratica, analizzarle in confronto a una configurazione passiva, e capire se l’algoritmo di ottimizzazione utilizzato fosse una tecnica efficiente per la taratura dei parametri. 
Nel corso del lavoro sono emersi però interessanti risultati e per questo si è voluto implementare un’ulteriore strategia di controllo, in modo da realizzare un’analisi più completa possibile.

Per la stesura della tesi si è fatto riferimento principalmente a Begnis et al., “An LMIbased approach for the control of semi-active magnetorheological suspensions”.
Nella realizzazione di questa tesi sono stati utilizzati in parallelo i software MATLAB e Simulink. Quest’ultimo necessario per la realizzazione dello schema a blocchi del modello quarter car e per le simulazioni del sistema dinamico, mentre con l’ausilio di MATLAB sono stati elaborati i parametri necessari per il sistema dinamico, è stato implementato l’algoritmo di ottimizzazione e sono stati raccolti i dati delle simulazioni.
