<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>MTBA: Matlab Toolbox for Biclustering Analysis - Contact</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="A Matlab toolbox for biclustering analysis">
    <meta name="author" content="Jayesh Kumar Gupta">

    <!--[if lt IE 9]>
        <script src="http://html5shim.googlecode.com/svn/trunk/html5.js"></script>
        <![endif]-->
    <link href="../css/bootstrap.min.css" rel="stylesheet">
    <link href="../css/bootstrap-responsive.min.css" rel="stylesheet">
    <link href="//netdna.bootstrapcdn.com/font-awesome/3.2.1/css/font-awesome.min.css" rel="stylesheet">
    <!-- add google analytics later -->
  </head>

  <body class="preview" id="top" data-spy="scroll" data-target=".subnav" data-offset="80">

    <!-- Navbar
         ================================================== -->
    <div class="navbar navbar-fixed-top">
      <div class="navbar-inner">
        <div class="container">
          <a class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </a>
          <a class="brand" href="http://iiirs.org/IIL/">Intelligent Informatics Lab</a>
          <div class="nav-collapse collapse" id="main-menu">
            <ul class="nav" id="main-menu-left">
              <li><a href="../people.html">People</a></li>
              <li><a href="../publications.html">Publications</a></li>
              <li><a href="../projects.html">Projects</a></li>
              <li class="dropdown" id="preview-menu">
                <a class="dropdown-toggle" data-toggle="dropdown" href="#">Download <b class="caret"></b></a>
                <ul class="dropdown-menu">
                  <li><a target="_blank" href="./">MTBA</a></li>
                  <li><a target="_blank" href="">Clustering</a></li>
                </ul>
              </li>
            </ul>

            <ul id="main-menu-right" class="nav pull-right">
              <li><a href="./contact.php">Contact Us</a></li>
              <li><form class="navbar-search pull-right" action="http://www.google.com/search" method="get">
                  <input type="hidden" name="sitesearch" value="http://iiirs.org/IIL/" />
                  <input name="q" type="text" class="search-query span2" placeholder="Search">
                </form>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </div>
    <div class="container">
      <br><br><br>
      <div class="content">
        <h1>Contact Us</h1>
        <div class="row">
          <div class="span8">
            <p>
        We are happy to receive feedback and bug reports or requests for more features, to discuss the toolbox in general as well as its documentation and to help you use it. You can also file issues at the <a href="https://bitbucket.org/rejuvyesh/matlab-toolbox-for-biclustering-analysis/issues?status=new&status=open">bitbucket repository</a>.
            </p>
            <p>
              We would also love to know how you use the toolbox while also hoping that you'll <a href="./#cite">cite us</a>.
            </p>
            <p>
              Contact us via the following form or mail us at <a href="mailto:iilatiitk@gmail.com">iilatiitk@gmail.com</a>.
            </p>
          </div>
        </div>
        <?php  
           
           // check for a successful form post  
           if (isset($_GET['s'])) echo "<div class=\"alert alert-success\">".$_GET['s']."</div>";  
           
           // check for a form error  
           elseif (isset($_GET['e'])) echo "<div class=\"alert alert-error\">".$_GET['e']."</div>";  
           
           ?>  
        <div class="row">
          <div class="span12">
            <form class="well span8 offset2" method="POST" action="contact-form-submission.php">
				<div class="row">
					<div class="span3">
						<label>First Name</label>
						<input type="text" name="first_name" class="span3" placeholder="Your First Name">
						<label>Last Name</label>
						<input type="text" name="last_name" class="span3" placeholder="Your Last Name">
						<label>Email Address</label>
						<div class="input-prepend">
							<span class="add-on"><i class="icon-envelope"></i></span><input type="text" id="inputIcon" name="contact_email" class="span2" style="width:233px" placeholder="Your email address">
						</div>
						<label>Subject
						  <select id="subject" name="subject" class="span3">
							  <option value="na" selected="">Choose One:</option>
							  <option value="general">General</option>
							  <option value="suggestions">Suggestions</option>
							  <option value="bugs">Bugs</option>
						  </select>
						</label>
					</div>
					<div class="span5">
						<label>Message</label>
						<textarea name="message" id="message" class="input-xlarge span5" rows="10"></textarea>
					</div>
				</div>
        <!-- <div class="form-actions">   -->
          <input type="hidden" name="save" value="contact">
				  <button type="submit" class="btn btn-primary pull-right">Send</button>
        <!-- </div>   -->
			      </form>
        </div></div>
      </div>
      
      <!-- Footer
           ================================================== -->
      <hr>

      <footer id="footer">
        <p class="pull-right"><a href="#top">Back to top</a></p>

        (c) Intelligent Informatics Lab<br/>
        Web fonts from <a href="http://www.google.com/webfonts">Google</a>. 
      </footer>

    </div><!-- /container -->



    <!-- Le javascript
         ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>
    <script src="../js/jquery.smooth-scroll.min.js"></script>
    <script src="../js/bootstrap.min.js"></script>
    <script src="../js/bootswatch.js"></script>


  </body>
</html>
