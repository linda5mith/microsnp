for dir in $(ls ./data); do
    echo "java -jar snpEff.jar build -genbank -v $dir" >> build_db_batch_cmds.sh;
done;
