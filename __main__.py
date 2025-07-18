import sys

from Sequential_Fish import viewer, pipeline, analysis, chromatic_abberrations, status
from Sequential_Fish._pipeline_scripts import PIPELINE_SCRIPTS
from Sequential_Fish.status import select_path_for_pipeline, select_path_for_analysis



def main():

    MODULES = ['viewer', 'pipeline', 'analysis', 'status', 'calibration']

    
    #CALL TO MODULES
    if len(sys.argv) < 2:
        print("Usage: python -m my_package <module> [args...]")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

    module = sys.argv[1]
    submodules = sys.argv[2:]

    if module == "viewer":

        if len(submodules) > 0 :
            print(f"No argument for viewer. Ignoring passed arguments : {submodules}")

        viewer.run()
        
    elif module == "pipeline":
        
        if '-p' in submodules : # for run in command line only : run path is passed after -p flag
            run_path_index = submodules.index('-p')
            run_path = submodules[run_path_index+1]
            submodules.pop(run_path_index)
            submodules.pop(run_path_index)
        else :
            run_path = select_path_for_pipeline() 

        if len(submodules) == 0 :
            pipeline.run(run_path) # This loads RUN_PATH from pipeline parameters and fix it for all scripts
        else :
            from Sequential_Fish.pipeline.runner import launch_script, script_folder # This loads RUN_PATH from pipeline parameters and fix it for all scripts
            if not all([script in PIPELINE_SCRIPTS for script in submodules]) :
                print(f"Unknown pipeline scripts. \nChoose from : {PIPELINE_SCRIPTS}")
            else :
                error_count = 0
                for script in submodules : 
                    if not script.endswith('.py') : script += ".py"
                    sucess = launch_script(script_folder + '/' + script, run_path=run_path)

                    if not sucess : error_count +=1
                
                print(f"all scripts ended with {error_count} errors. If any, error can be checked in run_log file.")
    
    elif module == "analysis" :

        if '-p' in submodules : # for run in command line only : run path is passed after -p flag
            run_path_index = submodules.index('-p')
            run_path = submodules[run_path_index+1]
            submodules.pop(run_path_index)
            submodules.pop(run_path_index)
        else :
            run_path = select_path_for_analysis() 
        
        if len(submodules) == 0 : 
            submodules = ['all']
            print("Starting all analysis modules")
        else :
            print("Starting selected analysis modules")
            
        analysis.run(run_path, *submodules)
        
        print("Done.")
    
    elif module == "status" :
        status.run()

    elif module == "calibration" :
        chromatic_abberrations.run()
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
        print("To skip path selection, use '-p' flag followed with fullpath to your directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()
