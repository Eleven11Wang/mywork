class brainspanwork 
    class sample_charactor(brainspanwork)
        function : filter the tissue that less then five 
        function : find which columns is which tissue (find_columns_tissue)
        function : seperate the matrix by age (seperate_by_time)
        seperate : the matrix by sex (seperate_by_sex)
    class expression_deal_with(sample_charactor)
        function : tissue expression matrix  return a work matrix(tissue_matrix_workon)
        function : 50% of the sample have this gene return a filted dataframe(filter_matrix_workon)
        function : check normalized or not ,box plot of data
        function : upper quantile transformation 
        function : log transform 
    
    (notdone)
    class differental_expression(sample_charactor)
        what i want to do first 
        analysis the differental expression of sex of whole matrix 
        analysis the differental expression of age of whole matrix 4 stage () 
        than filter by five  (rewrite it to a function) 
        tissue differental expression 
            tissue sex/age for 4 stage each stage/sex/spacial
    (workon)
    class basic_analysis_plot(sample_charactor)
        function : KDE plot of the matrix log2 expression distribution density
        function :average spearman coorelation of a region/area calculated for a period 
        def a expressed gene as a gene with at least one sample has expression value >6 and mean detection above in at least a period (DABG)
        euclidean distance (expression distance) 
            for sample all gene 495 495 2D plane plot 

        hierachial clustering ( sample tissue/age/sex )




    class expression_analysis(expression_deal_with)        
        function : differental expression age and sex 
        function : pca label by period/tissue (principle of pca how to plot that :w) 
        function : clustering 
        function : pearson correlation 
        what to deal with highexpression
    class network_analysis(expression_analysis):
        function : on class / paper research 

