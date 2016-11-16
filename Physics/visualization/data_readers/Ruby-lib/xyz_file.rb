require 'matrix'

module XYZ
  # Store the coordinates of each particle from the .xyz file into a matrix
  def getData(filename)
    mat = Matrix.empty(0,3)
    header = ""
    comment = ""

    filelines = File.open(filename).read.gsub(/\r\n?/, "\n")
    filelines.each_line.with_index do |line, i|
      if i==0
        header = line
        next
      elsif i==1
        comment = line
        next
      end

      numArray = line.scan(/[-+]?\d*\.?\d+/).map(&:to_f)
      mat = Matrix.rows(mat.to_a << numArray)
    end

    [header, comment, mat]
  end

  # Write data in DATAARRAY to FILENAME
  def putData(data, filename, *carriage)
    if carriage == [nil] || carriage == ["MAC"]
      carriage = "\n"
    elsif carriage == ["WIN"]
      carriage = "\r\n"
    end

    file = File.open(filename, 'w')
    file << "#{data[0]}#{data[1]}"
    data[2].row_vectors.each do |row|
      line = "H\t#{row[0]}\t#{row[1]}\t#{row[2]}#{carriage}"
      file << line
    end

    file.close()
  end
end
