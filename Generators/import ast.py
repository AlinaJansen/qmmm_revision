import ast

def extract_function_names(file_path):
    with open(file_path, 'r') as file:
        code = file.read()
        tree = ast.parse(code, filename=file_path)

    function_names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    
    for function_name in function_names:
        print(function_name)

# Passe den Pfad zu deiner Python-Datei an
file_path = '.\qmmm_revision\gmx2qmmm\pointcharges\sum_pcf_tm.py'
extract_function_names(file_path)