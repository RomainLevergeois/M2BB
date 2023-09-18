list_pckgs="argparse numpy matplotlib networkx Bio netgraph";

for element in $list_pckgs
do
    python3 -c "import $element";
    a=$(echo $?);
    if [ $a -ne 0 ]
    then
        echo "Installing $element :"
        pip install $element;
    else
        echo -e "$element already install\n";
    fi
done
echo "All the packages are install for TP1_languages_script.py";
